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

# 1. Create function to sample
# 2. Create function to run the simulation
# 3. Run simulation - fixed tmu (at prior), high N
# 4. Run simulation - fixed tmu (lower end), high N
# 5. Run simulation - fixed tmu (upper end), high N

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

    if (length(scaling) == 1) {
        scaling <- rep(scaling, n_sp)
    }

    # Generate occupancy probabilities based on SST
    q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2) * scaling[sid]

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
    dataset <- data.frame(sid = sid, sst = sst, y = y, prob_occ = occupancy, scaling = scaling[sid])

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



# 2. Create function to run the simulation -------
run_simulation <- function(tmus, tsds, N_pts, scaling, N_sims, stan_model_version = 5,
                           site_min = -Inf, site_max = Inf) {
    sim_results <- lapply(seq_len(N_sims), function(x) NULL)
    message("Using Stan model version ", stan_model_version, "\n")

    for (i in seq_len(N_sims)) {

        message(paste("Simulation", i, "of", N_sims))

        ds <- sim_dataset(
            N = length(tmus) * N_pts,
            n_sp = length(tmus),
            p = 0.5,
            tmu = tmus,
            tsd = rep(tsds, length(tmus)),
            scaling = scaling,
            site_min = site_min,
            site_max = site_max
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


# 3. Run simulation - fixed tmu (at prior), high N ------
# Number of simulations
N_sims <- 50
# Number of points/records
N_pts <- 500
# Means of sites -> prior of STAN model is 20
tmus <- rep(20, 5)
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- seq(0.2, 1, length.out = 5)


sim1_v4 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 4
)
sim1_v5 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 5
)

# Pull 1 simulation as example
sim1_v4 |> filter(simulation_id == 1)
sim1_v5 |> filter(simulation_id == 1)

plot_histogram(sim1_v4, "Average value - with scaling - Model 4")
plot_histogram(sim1_v5, "Average value - with scaling - Model 5")




# 4. Run simulation - fixed tmu (lower end), high N ------
# Number of simulations
N_sims <- 50
# Number of points/records
N_pts <- 500
# Means of sites -> prior of STAN model is 20
tmus <- rep(15, 5)
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- seq(0.2, 1, length.out = 5)


sim2_v4 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 4
)
sim2_v5 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 5
)

# Pull 1 simulation as example
sim2_v4 |> filter(simulation_id == 1)
sim2_v5 |> filter(simulation_id == 1)

plot_histogram(sim2_v4, "Average value - with scaling - lower end - Model 4")
plot_histogram(sim2_v5, "Average value - with scaling - lower end - Model 5")




# 5. Run simulation - fixed tmu (upper end), high N ------
# Number of simulations
N_sims <- 50
# Number of points/records
N_pts <- 500
# Means of sites -> prior of STAN model is 20
tmus <- rep(28, 5)
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- seq(0.2, 1, length.out = 5)


sim3_v4 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 4
)
sim3_v5 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 5
)

# Pull 1 simulation as example
sim3_v4 |> filter(simulation_id == 1)
sim3_v5 |> filter(simulation_id == 1)

plot_histogram(sim3_v4, "Average value - with scaling - upper end - Model 4")
plot_histogram(sim3_v5, "Average value - with scaling - upper end - Model 5")



# 6. Run simulation - varying tmu, high N ------
# Number of simulations
N_sims <- 50
# Number of points/records
N_pts <- 500
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- sample(seq(0.2, 1, length.out = 5), 15, replace = T)


sim4_v4 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 4
)
sim4_v5 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 5
)

# Pull 1 simulation as example
sim4_v4 |> filter(simulation_id == 1)
sim4_v5 |> filter(simulation_id == 1)

plot_histogram(sim4_v4, "Average value - with scaling - Model 4")
plot_histogram(sim4_v5, "Average value - with scaling - Model 5")



# 7. Run simulation - varying tmu, different N points ------
# Number of simulations
N_sims <- 1
# Number of points/records
N_pts <- c(seq(500, 100, by = -100), 20)
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- sample(seq(0.2, 1, length.out = 5), 15, replace = T)


sim5_v4 <- lapply(N_pts, \(x) {
    message("Running with N pts = ", x)
    sim <- run_simulation(
        tmus = tmus, tsds = tsds, N_pts = x, N_sims = N_sims,
        scaling = rare_scaling, stan_model_version = 4
    )
    sim$N_pts <- x
    sim
})
sim5_v5 <- lapply(N_pts, \(x) {
    message("Running with N pts = ", x)
    sim <- run_simulation(
        tmus = tmus, tsds = tsds, N_pts = x, N_sims = N_sims,
        scaling = rare_scaling, stan_model_version = 5
    )
    sim$N_pts <- x
    sim
})

sim5_v4 <- do.call("rbind", sim5_v4)
sim5_v5 <- do.call("rbind", sim5_v5)

sim5_v4 |> 
    filter(what == "tmu") |>
    ggplot(aes(x = expected, y = mean)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~ N_pts) + theme_light()

sim5_v5 |> 
    filter(what == "tmu") |>
    ggplot(aes(x = expected, y = mean)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(~ N_pts) + theme_light()



# 8. Run simulation - varying tmu, high N, truncated at high point ------
# Number of simulations
N_sims <- 50
# Number of points/records
N_pts <- 500
# Means of sites -> prior of STAN model is 20
tmus <- 16:30 # Different from the others, to force more the truncation on the last
# SD (fixed across sites)
tsds <- 2
# Set scaling (just for second simulation)
rare_scaling <- sample(seq(0.2, 1, length.out = 5), 15, replace = T)
# Maximum of site for truncation
site_max <- 31

sim6_v4 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 4, site_max = site_max
)
sim6_v5 <- run_simulation(
    tmus = tmus, tsds = tsds, N_pts = N_pts, N_sims = N_sims,
    scaling = rare_scaling, stan_model_version = 5, site_max = site_max
)

# Pull 1 simulation as example
sim6_v4 |> filter(simulation_id == 1)
sim6_v5 |> filter(simulation_id == 1)

plot_histogram(sim6_v4, "Average value - with scaling - Model 4")
plot_histogram(sim6_v5, "Average value - with scaling - Model 5")



# 9. Real data
library(terra)
species <- c("boru", "hara", "mybo", "lujo", "spam", "scze", "acch")
sst <- rast("data-raw/thetao_baseline_depthsurf_mean.tif")

species_data <- lapply(species, function(x){
    dat <- read.csv(file.path("data", paste0(x, "_pa.csv")))
    colnames(dat)[3] <- "presence"
    dat$species <- x
    sst_d <- extract(sst, dat[,1:2])
    dat$sst <- sst_d$mean
    dat <- dat[!is.na(dat$sst),]
    dat
})

species_data <- do.call("rbind", species_data)

expected <- species_data |> 
    filter(presence == 1) |> 
    group_by(species) |> 
    summarise(expected = mean(sst), expected_sd = sd(sst))

species_data <- left_join(species_data, expected)

dat <- list(
            N = nrow(species_data),
            N_spp = length(unique(species_data$species)),
            sid = as.integer(as.factor(species_data$species)),
            sst = species_data$sst,
            y = species_data$presence,
            expected = species_data$expected,
            expected_sd = species_data$expected_sd
        )

m5 <- cstan(file = paste0("codes/model5.stan"), data = dat, rstan_out = FALSE)
prec_res_m5 <- precis(m5, 2)

m4 <- cstan(file = paste0("codes/model4.stan"), data = dat, rstan_out = FALSE)
prec_res_m4 <- precis(m4, 2)

plot_data <- data.frame(
    species_id = c(1:7, 1:7),
    predicted = c(prec_res_m4$mean[1:7], prec_res_m5$mean[1:7]),
    predicted_low = c(prec_res_m4$`5.5%`[1:7], prec_res_m5$`5.5%`[1:7]),
    predicted_high = c(prec_res_m4$`94.5%`[1:7], prec_res_m5$`94.5%`[1:7]),
    model = c(rep("M4", 7), rep("M5", 7))
)

expected2 <- expected
expected2$species_id <- as.integer(as.factor(expected2$species))

plot_data <- left_join(plot_data, expected2)

ggplot(plot_data) +
    geom_pointrange(aes(x = expected, y = predicted, ymin = predicted_low, ymax = predicted_high, color = species), linewidth = 1) +
    geom_pointrange(aes(x = expected, y = predicted, xmin = expected - expected_sd, xmax = expected + expected_sd, color = species),
    linetype = 2) +
    theme_light() +
    scale_x_continuous(breaks = seq(14, 32, by = 2), limits = c(14,32)) +
    scale_y_continuous(breaks = seq(14, 32, by = 2), limits = c(14,32)) +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "black")) +
    facet_wrap(~model)
ggsave("figures/true_fish.png", width = 12, height = 7)


wrld <- rnaturalearth::ne_countries(returnclass = "sf")

species_data_sf <- sf::st_as_sf(species_data, coords = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
sf::sf_use_s2(FALSE)
wrld <- sf::st_crop(wrld, species_data_sf)
species_data_sf$presence <- as.factor(species_data_sf$presence)
ggplot() +
    geom_sf(data = wrld, fill = "grey70", color = "grey60") +
    geom_sf(data = species_data_sf, aes(color = presence), alpha = .5) +
    scale_color_manual(values = c("#1132b7", "#f57d05")) +
    theme_light() +
    facet_wrap(~ species)
ggsave("figures/true_fish_maps.png", width = 16, height = 11)

save.image(file = "to-continue.RData")
load("to-continue.RData")
