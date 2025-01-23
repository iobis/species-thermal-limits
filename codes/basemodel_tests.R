# Base model simulation
set.seed(2023)
library(rethinking)
library(dplyr)
library(ggplot2)

# Create a helper function to simulate datasets
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

# Create other helper functions
plot_histogram <- function(sim_results_agg, mode = "simple") {
    if (mode == "simple") {
        sim_results_agg |>
            filter(what == "tmu") |>
            ggplot() +
                geom_histogram(aes(x = mean)) +
                geom_vline(aes(xintercept = base_val), color = "#cb8917") +
                facet_wrap(~species) +
                theme_light() +
                theme(panel.grid = element_blank()) +
                ggtitle("Average value - 100 simulations")
    } else {
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
                ggtitle("Average value - 100 simulations")
    }
}

plot_single_sim <- function(sim_results_agg, sim = 1, which_data = "tmu") {
    sim_1 <- sim_results_agg |>
    filter(what == which_data, simulation_id == sim)
    with(sim_1,
        plot(x = expected, y = mean, col = "#21568f", pch = 20, cex = 2,
        ylim = range(c(expected, mean))))
    with(sim_1,
        points(x = base_val, y = mean, col = "#b5c61c71", pch = 20, cex = 2))
    abline(with(sim_1, lm(mean ~ expected)))
    legend("topleft", legend = c("Mean of sample", "Initial mean"),
       col = c("#21568f", "#b5c61c71"), pch = 20)
}

pull_stats <- function(sim_datasets, sim = 1, which_stat = "count", format = "wider") {
    ds <- sim_datasets[[i]]$dataset

    if (which_stat == "count") {
        result <- ds |> group_by(sid, y) |> count()
    } else {
        result <- ds |> group_by(sid, y) |> 
            summarise(mean = mean(sst), sd = sd(sst), mean_prob = mean(prob_occ))
    }
    if (format == "wider") {
        result <- result |> 
            tidyr::pivot_wider(names_from = y, values_from = !all_of(c("sid", "y")))
    }
    return(result)
}


#### Run simulations
# SIMULATION 1
# Normal conditions
# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2

sim_results <- lapply(seq_len(N_sims), function(x) NULL)
sim_datasets <- sim_results

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = 0.5,
        tmu = tmus,
        tsd = rep(tsds, length(tmus))
    )
    sim_datasets[[i]] <- ds

    dat <- list(
        N = nrow(ds$dataset),
        N_spp = length(unique(ds$dataset$sid)),
        sid = as.integer(as.factor(ds$dataset$sid)),
        sst = ds$dataset$sst,
        y = ds$dataset$y
    )
    m <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

    prec_res <- precis(m, 2)
    prec_res$expected <- c(
        unname(tapply(ds$dataset$sst, ds$dataset$sid, mean)),
        unname(tapply(ds$dataset$sst, ds$dataset$sid, sd)),
        0
    )
    prec_res$base_val <- c(tmus, rep(tsds, length(tmus)), 0)
    prec_res$species <- c(seq_along(tmus), seq_along(tmus), 0)
    prec_res$what <- c(rep("tmu", length(tmus)), rep("tsd", length(tmus)), "p")
    prec_res$simulation_id <- i

    sim_results[[i]] <- prec_res
}

sim_results_agg <- do.call("rbind", sim_results)

plot_histogram(sim_results_agg, mode = "simple")
plot_histogram(sim_results_agg, mode = "complete")

# Separate one simulation
plot_single_sim(sim_results_agg, sim = 1)
pull_stats(sim_datasets, sim = 1, which_stat = "count")
pull_stats(sim_datasets, sim = 1, which_stat = "means")


# SIMULATION 2
# Truncated at upper bound
# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2
# Define a maximum for site (truncated)
max_site <- 31

sim2_results <- lapply(seq_len(N_sims), function(x) NULL)
sim2_datasets <- sim2_results

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = 0.5,
        tmu = tmus,
        tsd = rep(tsds, length(tmus)),
        site_max = max_site
    )
    sim2_datasets[[i]] <- ds

    dat <- list(
        N = nrow(ds$dataset),
        N_spp = length(unique(ds$dataset$sid)),
        sid = as.integer(as.factor(ds$dataset$sid)),
        sst = ds$dataset$sst,
        y = ds$dataset$y
    )
    m <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

    prec_res <- precis(m, 2)
    prec_res$expected <- c(
        unname(tapply(ds$dataset$sst, ds$dataset$sid, mean)),
        unname(tapply(ds$dataset$sst, ds$dataset$sid, sd)),
        0
    )
    prec_res$base_val <- c(tmus, rep(tsds, length(tmus)), 0)
    prec_res$species <- c(seq_along(tmus), seq_along(tmus), 0)
    prec_res$what <- c(rep("tmu", length(tmus)), rep("tsd", length(tmus)), "p")
    prec_res$simulation_id <- i

    sim2_results[[i]] <- prec_res
}

sim2_results_agg <- do.call("rbind", sim2_results)

plot_histogram(sim2_results_agg, mode = "simple")
plot_histogram(sim2_results_agg, mode = "complete")

# Separate one simulation
plot_single_sim(sim2_results_agg, sim = 1)
pull_stats(sim2_datasets, sim = 1, which_stat = "count")
pull_stats(sim2_datasets, sim = 1, which_stat = "means")





# SIMULATION 3
# Truncated at upper bound and with a lower detectability for absence
# i.e., chances of wrongly identifying an absence are higher
# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2
# Detectability (only for absence)
detect <- 0.1
# Define a maximum for site (truncated)
max_site <- 31

sim3_results <- lapply(seq_len(N_sims), function(x) NULL)
sim3_datasets <- sim3_results

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = detect,
        tmu = tmus,
        tsd = rep(tsds, length(tmus)),
        site_max = max_site
    )
    sim3_datasets[[i]] <- ds

    dat <- list(
        N = nrow(ds$dataset),
        N_spp = length(unique(ds$dataset$sid)),
        sid = as.integer(as.factor(ds$dataset$sid)),
        sst = ds$dataset$sst,
        y = ds$dataset$y
    )
    m <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

    prec_res <- precis(m, 2)
    prec_res$expected <- c(
        unname(tapply(ds$dataset$sst, ds$dataset$sid, mean)),
        unname(tapply(ds$dataset$sst, ds$dataset$sid, sd)),
        0
    )
    prec_res$base_val <- c(tmus, rep(tsds, length(tmus)), 0)
    prec_res$species <- c(seq_along(tmus), seq_along(tmus), 0)
    prec_res$what <- c(rep("tmu", length(tmus)), rep("tsd", length(tmus)), "p")
    prec_res$simulation_id <- i

    sim3_results[[i]] <- prec_res
}

sim3_results_agg <- do.call("rbind", sim3_results)

plot_histogram(sim3_results_agg, mode = "simple")
plot_histogram(sim3_results_agg, mode = "complete")

# Separate one simulation
plot_single_sim(sim3_results_agg, sim = 1)
pull_stats(sim3_datasets, sim = 1, which_stat = "count")
pull_stats(sim3_datasets, sim = 1, which_stat = "means")




# SIMULATION 4
# In this case we apply a scaling constant to reduce the maximum probability of
# encountering the species to a very low value. The species is rare and thus
# there will be many absences even in places that are "suitable"
# In terms of the gaussian probability/suitability, this is the behavior:
prob_curve <- function(tmu, tsd, sst = 25, term = 1) {
    exp(-0.5 * ((sst - tmu) / tsd)^2) * term
}
par(mfrow = c(1,1))
sst_vals <- seq(20, 30, length.out = 100)
normal_curve <- prob_curve(sst_vals, 2)
plot(normal_curve ~ sst_vals, ylab = "Probability", xlab = "SST")
rare_curve <- prob_curve(sst_vals, 2, term = 0.5)
points(rare_curve ~ sst_vals, col = "blue")
legend("topright", legend = c("Normal", "Rescaled"),
       col = c("black", "blue"), pch = 1)

# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2
# Detectability (only for absence)
detect <- 0.5
# Scaling factor
scaling <- 0.2


sim4_results <- lapply(seq_len(N_sims), function(x) NULL)
sim4_datasets <- sim4_results

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = detect,
        tmu = tmus,
        tsd = rep(tsds, length(tmus)),
        site_max = 31,
        scaling = scaling
    )
    sim4_datasets[[i]] <- ds

    dat <- list(
        N = nrow(ds$dataset),
        N_spp = length(unique(ds$dataset$sid)),
        sid = as.integer(as.factor(ds$dataset$sid)),
        sst = ds$dataset$sst,
        y = ds$dataset$y
    )
    m <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

    prec_res <- precis(m, 2)
    prec_res$expected <- c(
        unname(tapply(ds$dataset$sst, ds$dataset$sid, mean)),
        unname(tapply(ds$dataset$sst, ds$dataset$sid, sd)),
        0
    )
    prec_res$base_val <- c(tmus, rep(tsds, length(tmus)), 0)
    prec_res$species <- c(seq_along(tmus), seq_along(tmus), 0)
    prec_res$what <- c(rep("tmu", length(tmus)), rep("tsd", length(tmus)), "p")
    prec_res$simulation_id <- i

    sim4_results[[i]] <- prec_res
}

sim4_results_agg <- do.call("rbind", sim4_results)

plot_histogram(sim4_results_agg, mode = "simple")
plot_histogram(sim4_results_agg, mode = "complete")

# Separate one simulation
plot_single_sim(sim4_results_agg, sim = 1)
pull_stats(sim4_datasets, sim = 1, which_stat = "count")
pull_stats(sim4_datasets, sim = 1, which_stat = "means")
