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
    site_max = Inf                              # Maximum of site (for truncation)
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
    q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2)

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



#### Run simulations
# SIMULATION 1
# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2

sim_results <- lapply(seq_len(N_sims), function(x) NULL)

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = 0.5,
        tmu = tmus,
        tsd = rep(tsds, length(tmus))
    )
    #data.frame(true_mean = tmus, sim_means = tapply(ds$dataset$sst, ds$dataset$sid, mean))

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

sps <- seq_along(tmus)

sim_results_agg |>
    filter(what == "tmu") |>
    ggplot() +
        geom_histogram(aes(x = mean)) +
        geom_vline(aes(xintercept = base_val), color = "#cb8917") +
        facet_wrap(~species) +
        theme_light() +
        theme(panel.grid = element_blank()) +
        ggtitle("Average value - 100 simulations")

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




# SIMULATION 2
# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2

sim2_results <- lapply(seq_len(N_sims), function(x) NULL)

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = 0.5,
        tmu = tmus,
        tsd = rep(tsds, length(tmus)),
        site_max = 31
    )
    #data.frame(true_mean = tmus, sim_means = tapply(ds$dataset$sst, ds$dataset$sid, mean))

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

sim2_results_agg |>
    filter(what == "tmu") |>
    ggplot() +
        geom_histogram(aes(x = mean)) +
        geom_vline(aes(xintercept = base_val), color = "#cb8917") +
        facet_wrap(~species) +
        theme_light() +
        theme(panel.grid = element_blank()) +
        ggtitle("Average value - 100 simulations")

sim2_results_agg |>
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

# Separate one simulation
sim_1 <- sim2_results_agg |>
    filter(what == "tmu", simulation_id == 1)
with(sim_1,
     plot(x = expected, y = mean, col = "#21568f", pch = 20, cex = 2))
with(sim_1,
     points(x = base_val, y = mean, col = "#b5c61c71", pch = 20, cex = 2))
abline(with(sim_1, lm(mean ~ expected)))




# SIMULATION 3
# Number of simulations
N_sims <- 100
# Means of sites -> prior of STAN model is 20
tmus <- 15:29
# SD (fixed across sites)
tsds <- 2
# Detectability (only for absence)
detect <- 0.1


sim3_results <- lapply(seq_len(N_sims), function(x) NULL)

for (i in seq_len(N_sims)) {

    message(paste("Simulation", i, "of", N_sims))

    ds <- sim_dataset(
        N = length(tmus) * 100,
        n_sp = length(tmus),
        p = detect,
        tmu = tmus,
        tsd = rep(tsds, length(tmus)),
        site_max = 31
    )
    #data.frame(true_mean = tmus, sim_means = tapply(ds$dataset$sst, ds$dataset$sid, mean))

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

sim3_results_agg |>
    filter(what == "tmu") |>
    ggplot() +
        geom_histogram(aes(x = mean)) +
        geom_vline(aes(xintercept = base_val), color = "#cb8917") +
        facet_wrap(~species) +
        theme_light() +
        theme(panel.grid = element_blank()) +
        ggtitle("Average value - 100 simulations")

sim3_results_agg |>
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

# Separate one simulation
sim_1 <- sim3_results_agg |>
    filter(what == "tmu", simulation_id == 1)
with(sim_1,
     plot(x = expected, y = mean, col = "#21568f", pch = 20, cex = 2))
with(sim_1,
     points(x = base_val, y = mean, col = "#b5c61c71", pch = 20, cex = 2))
abline(with(sim_1, lm(mean ~ expected)))
