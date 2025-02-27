# Base model simulation (minimal example)

# Load packages
set.seed(2023)
library(rethinking)
library(dplyr)
library(ggplot2)
# For truncated normal
if (!require(truncnorm, quietly = T)) {
    install.packages("truncnorm")
}

# 1. Create function to sample (line 19)
# 2. Set the parameters for the simulation (line 102)
# 3. Create function to run the simulation (line 115)
# 4. Run simulations with no scaling (line 170)
# 5. Run simulation with scaling, simulating low probability of presence (line 184)

# 1. Create function to sample ------
# We create a function to generate the simulated datasets
# Our simulated model follows exactly the same process from the model
# i.e. we generate occupancy probabilities using:
# q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2) * scaling
# being `sst` the SST vector, `tmu` the mean sst for species `sid` and `tsd` the 
# SD for species `sid`. We also add a scaling factor, that is 1 for the normal
# conditions (that is, no modification) or any other value < 1, what will
# cause the maximum probability to be lower than 1 (squeeze the gaussian curve)
# This simulates cases that a species is rare
#
# With the occupancy probabilities, we then sample presence/absence.
# We assume that for presence (1) the probability is exactly the occupancy probability
# while for absence there is the chance that is a false absence, given by:
# prob_absent <- (1 - prob_occupied) + prob_occupied * (1 - p)
sim_dataset <- function(
    N = 500,                                    # Number of surveys
    n_sp = 5,                                   # Number of species
    p = rbeta(1, 2, 2),                         # Detection probability - fixed value
    tmu = round(rnorm(n_sp, 20, 3), 1),         # Mean SST per species
    tsd = round(abs(rnorm(n_sp, 5, 1)), 1),     # SD of SST per species
    site_min = -Inf,                            # Minimum of site (for truncation)
    site_max = Inf,                             # Maximum of site (for truncation)
    scaling = 1                                 # Scaling to modify maximum probability
) {
    # Simulate surveys
    sid <- rep(seq_len(n_sp), each = N/n_sp)    # Species ID
    if (site_min != -Inf || site_max != Inf) {
        message("Unsing truncated normal")
        sst <- truncnorm::rtruncnorm(
            n = N, a = site_min, b = site_max, mean = tmu[sid], sd = tsd[sid]
        ) # Observed SST across surveys (truncated)
    } else {
        sst <- rnorm(N, tmu[sid], tsd[sid])   # Observed SST across surveys
    }

    # Generate occupancy probabilities based on SST
    q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2) * scaling

    # Simulate observations
    y <- numeric(N)
    occupancy  <- numeric(N)
    for (i in 1:N) {
        prob_occupied <- q[i] # Probability of occupancy given suitability
        occupancy[i]  <- prob_occupied
        # Detected if occupied
        prob_detected <- prob_occupied #presence is always true, otherwise use next line
        # prob_detected <- prob_occupied * p
        # Absence or undetected
        prob_absent <- (1 - prob_occupied) + prob_occupied * (1 - p)
        y[i] <- sample(c(1, 0), size = 1, prob = c(prob_detected, prob_absent))
    }

    # Create dataset
    dataset <- data.frame(sid = sid, sst = sst, y = y, prob_occ = occupancy)

    list(
        p = p,
        tmu = tmu,
        tsd = tsd,
        dataset = dataset
    )
}

# We also create one more function to make the histogram plots
plot_histogram <- function(sim_results_agg, title) {
    sim_results_agg |>
            filter(what == "tmu") |>
            ggplot() +
                geom_histogram(aes(x = `5.5%`), fill = "#0044b9", alpha = 0.5) +
                geom_histogram(aes(x = `94.5%`), fill = "#a5004d", alpha = 0.5) +
                geom_histogram(aes(x = mean), fill = "grey50") +
                geom_vline(aes(xintercept = base_val), color = "#000000", linewidth = 1.2) +
                geom_vline(aes(xintercept = base_val - (1.96 * 2)), color = "#7c7c7c", linewidth = 1) +
                geom_vline(aes(xintercept = base_val + (1.96 * 2)), color = "#7c7c7c", linewidth = 1) +
                facet_wrap(~species) +
                theme_light() +
                theme(panel.grid = element_blank()) +
                ggtitle(title)
}




# 2. Set the parameters for the simulation -----
# Number of simulations
N_sims <- 20
# Means of sites -> prior of STAN model is 20
tmus <- 15:29 # 15 species (each with an optimum SST)
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- 0.2



# 3. Create function to run the simulation -------
run_simulation <- function(scaling = 1, stan_model_version = 5) {
    sim_results <- lapply(seq_len(N_sims), function(x) NULL)
    message("Using Stan model version ", stan_model_version, "\n")

    for (i in seq_len(N_sims)) {

        message(paste("Simulation", i, "of", N_sims))

        ds <- sim_dataset(
            N = length(tmus) * 100,
            n_sp = length(tmus),
            p = 0.5,
            tmu = tmus,
            tsd = rep(tsds, length(tmus))
        )

        dat <- list(
            N = nrow(ds$dataset),
            N_spp = length(unique(ds$dataset$sid)),
            sid = as.integer(as.factor(ds$dataset$sid)),
            sst = ds$dataset$sst,
            y = ds$dataset$y
        )
        m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = dat, rstan_out = FALSE)

        prec_res <- precis(m, 2)
        if (stan_model_version != 4) {
            to_fill <- length(tmus)
            what <- rep("tomax", length(tmus))
            to_repeat <- 3
        } else {
            to_fill <- 0
            to_repeat <- 2
            what <- NULL
        }
        prec_res$expected <- c(
            unname(tapply(ds$dataset$sst, ds$dataset$sid, mean)),
            unname(tapply(ds$dataset$sst, ds$dataset$sid, sd)),
            rep(NA, to_fill + 1)
        )
        prec_res$base_val <- c(tmus, rep(tsds, length(tmus)), rep(NA, to_fill + 1))
        prec_res$species <- c(rep(seq_along(tmus), to_repeat), NA)
        prec_res$what <- c(rep("tmu", length(tmus)), rep("tsd", length(tmus)), what, "p")
        prec_res$simulation_id <- i

        sim_results[[i]] <- prec_res
    }

    do.call("rbind", sim_results)
}




# 4. Run simulations with no scaling (i.e. normal conditions) ------
sim_results_agg_v4 <- run_simulation(stan_model_version = 4)
sim_results_agg_v5 <- run_simulation(stan_model_version = 5)

# Pull 1 simulation as example
sim_results_agg_v4 |> filter(simulation_id == 1)
sim_results_agg_v5 |> filter(simulation_id == 1)

plot_histogram(sim_results_agg_v4, "Average value - no scaling - Model 4")
plot_histogram(sim_results_agg_v5, "Average value - no scaling - Model 5")
# Plots shows the distribution of estimates, for each species, across simulations
# The blue histogram depicts the 5.5% value, the red the 94.5%, while the grey represents the mean.
# Black line is the original mean for the species, and grey lines mean+SD

# 5. Run simulation with scaling, simulating low probability of presence ------
sim_results_agg_rs_v4 <- run_simulation(scaling = rare_scaling, stan_model_version = 4)
sim_results_agg_rs_v5 <- run_simulation(scaling = rare_scaling, stan_model_version = 5)

# Pull 1 simulation as example
sim_results_agg_rs_v4 |> filter(simulation_id == 1)
sim_results_agg_rs_v5 |> filter(simulation_id == 1)

plot_histogram(sim_results_agg_rs_v4, "Average value - with scaling - Model 4")
plot_histogram(sim_results_agg_rs_v5, "Average value - with scaling - Model 5")
