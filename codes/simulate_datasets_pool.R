sim_dataset_pool <- function(
    N = 500,                                    # Number of surveys
    n_sp = 5,                                   # Number of species
    p = rbeta(1, 2, 2),                         # Detection probability
    mu_tmu = 20,                                # Global mean for tmu
    sigma_tmu = 3,                              # Global sd for tmu
    mu_tsd = 5,                                 # Global mean for tsd
    sigma_tsd = 1,                              # Global sd for tsd 
    tomax = rbeta(n_sp, 5, 1),                  # Max occupancy prob per species - you can also pass a fixed value
    site_min = -Inf,                            # Minimum of site (for truncation)
    site_max = Inf                              # Maximum of site (for truncation)
) {

    tmu <- round(rnorm(n_sp, mu_tmu, sigma_tmu), 1)
    tsd <- round(abs(rnorm(n_sp, mu_tsd, sigma_tsd)), 1)

    # Species IDs
    sid <- rep(seq_len(n_sp), each = N / n_sp)

    # Simulate SST for each survey
    if (site_min != -Inf || site_max != Inf) {
        sst <- truncnorm::rtruncnorm(n = N, a = site_min, b = site_max, mean = tmu[sid], sd = tsd[sid])
    } else {
        sst <- rnorm(N, mean = tmu[sid], sd = tsd[sid])
    }

    if (length(tomax) == 1) {
        tomax <- rep(tomax, n_sp)
    }

    # Calculate occupancy probability (Gaussian suitability * tomax)
    q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2) * tomax[sid]

    # Simulate detections
    y <- numeric(N)
    occupancy <- numeric(N)
    for (i in 1:N) {
        prob_occupied <- q[i] # Probability of occupancy given suitability
        occupancy[i] <- prob_occupied
        # Detected if occupied
        # prob_detected <- prob_occupied #presence is always true, otherwise use next line
        # prob_detected <- prob_occupied * p
        # Absence or undetected
        prob_absent <- (1 - prob_occupied) + prob_occupied * (1 - p)
        y[i] <- sample(c(1, 0), size = 1, prob = c(prob_occupied * p, prob_absent))
    }

    # Create dataset
    dataset <- data.frame(sid = sid, sst = sst, y = y, prob_occ = occupancy)
    list(
        dataset = dataset,
        p = p,
        tmu = tmu,
        tsd = tsd,
        tomax = tomax,
        mu_tmu = mu_tmu,
        sigma_tmu = sigma_tmu,
        mu_tsd = mu_tsd,
        sigma_tsd = sigma_tsd
    )
}

# Example usage
set.seed(2025)
library(rethinking)
library(dplyr)
stan_model_version <- 6

sim_data <- sim_dataset_pool(
    N = 2000,
    n_sp = 20,
    p = 0.5,
    mu_tmu = 25,
    sigma_tmu = 3,
    mu_tsd = 5,
    sigma_tsd = 1,
    site_max =  31
)

dat <- list(
    N = nrow(sim_data$dataset),
    N_spp = length(unique(sim_data$dataset$sid)),
    sid = as.integer(as.factor(sim_data$dataset$sid)),
    sst = sim_data$dataset$sst,
    y = sim_data$dataset$y
)

m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = dat, rstan_out = FALSE)

prec_res <- precis(m, 2)
prec_res$expected <- c(sim_data$tmu, sim_data$tsd, sim_data$tomax,
                       sim_data$p, sim_data$mu_tmu, sim_data$sigma_tmu,
                       sim_data$mu_tsd, sim_data$sigma_tsd)
prec_res$delta <- prec_res$mean - prec_res$expected
prec_res

# Now try to vary number of points per species
sim_data_b <- sim_data
for (i in unique(sim_data_b$dataset$sid)) {
    isamp <- (\(){
        v <- as.integer(rnorm(1, 70, 20))
        ifelse(v < 10, 10, ifelse(v > 100, 100, v))
    })()
    sim_data_b$dataset$to_remove[
        sim_data_b$dataset$sid == i
    ] <- sample(c(rep(TRUE, isamp), rep(FALSE, 100-isamp)))
}
sim_data_b$dataset <- sim_data_b$dataset[sim_data_b$dataset$to_remove,]
sim_data_b$n_sampled <- sim_data_b$dataset |>
    group_by(sid) |>
    summarise(n_sampled = n()) |> pull()

dat <- list(
    N = nrow(sim_data_b$dataset),
    N_spp = length(unique(sim_data_b$dataset$sid)),
    sid = as.integer(as.factor(sim_data_b$dataset$sid)),
    sst = sim_data_b$dataset$sst,
    y = sim_data_b$dataset$y
)

m_b <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = dat, rstan_out = FALSE)

prec_res_b <- precis(m_b, 2)
prec_res_b$expected <- c(sim_data_b$tmu, sim_data_b$tsd, sim_data_b$tomax,
                       sim_data_b$p, sim_data_b$mu_tmu, sim_data_b$sigma_tmu,
                       sim_data_b$mu_tsd, sim_data_b$sigma_tsd)
prec_res_b$delta <- prec_res_b$mean - prec_res_b$expected
prec_res_b
plot(prec_res_b$mean[grep("tmu\\[", rownames(prec_res_b))],
     prec_res_b$expected[grep("tmu\\[", rownames(prec_res_b))],
     col = "#1c57be6b", pch = 19,
     xlab = "Estimated tmu", ylab = "Expected tmu",
     main = "Estimated vs Expected tmu",
     xlim = c(range(c(prec_res_b$mean[grep("tmu\\[", rownames(prec_res_b))],
                      prec_res_b$expected[grep("tmu\\[", rownames(prec_res_b))]))),
     ylim = c(range(c(prec_res_b$mean[grep("tmu\\[", rownames(prec_res_b))],
                      prec_res_b$expected[grep("tmu\\[", rownames(prec_res_b))]))))
# Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
# Chain 1 Exception: normal_lpdf: Scale parameter is 0, but must be positive! (in '/var/folders/jd/hdd9v5xd2zs_296j2ksd9px80000gp/T/RtmpfPuBo6/model-1255364556076.stan', line 32, column 4 to column 36)
# Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
# Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
# Chain 1 
#Warning: 2 of 500 (0.0%) transitions ended with a divergence.
#See https://mc-stan.org/misc/warnings for details.

# Now do test with multiple simulations to see how consistently is getting the expected values
n_simulations <- 30
results_list <- vector("list", n_simulations)
for (s in seq_len(n_simulations)) {
    message("\nSimulation ", s, "\n\n")

    sim_data <- sim_dataset_pool(
        N = 2000,
        n_sp = 20,
        p = 0.5,
        mu_tmu = 20,
        sigma_tmu = 3,
        mu_tsd = 5,
        sigma_tsd = 1,
        site_max =  31
    )

    dat <- list(
        N = nrow(sim_data$dataset),
        N_spp = length(unique(sim_data$dataset$sid)),
        sid = as.integer(as.factor(sim_data$dataset$sid)),
        sst = sim_data$dataset$sst,
        y = sim_data$dataset$y
    )

    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = dat, rstan_out = FALSE)

    prec_res <- precis(m, 2)
    prec_res$expected <- c(sim_data$tmu, sim_data$tsd, sim_data$tomax,
                        sim_data$p, sim_data$mu_tmu, sim_data$sigma_tmu,
                        sim_data$mu_tsd, sim_data$sigma_tsd)
    prec_res$delta <- prec_res$mean - prec_res$expected
    results_list[[s]] <- prec_res
}

results_list <- do.call(rbind, results_list)

tmus <- results_list$mean[grep("tmu\\[", rownames(results_list))]
tmus_e <- results_list$expected[grep("tmu\\[", rownames(results_list))]

plot(tmus, tmus_e, col = "#1c57be6b", pch = 19,
     xlab = "Estimated tmu", ylab = "Expected tmu",
     main = "Estimated vs Expected tmu",
     xlim = c(range(c(tmus, tmus_e))),
     ylim = c(range(c(tmus, tmus_e))))

tsds <- results_list$mean[grep("tsd\\[", rownames(results_list))]
tsds_e <- results_list$expected[grep("tsd\\[", rownames(results_list))]

plot(tsds, tsds_e, col = "#1c57be6b", pch = 19,
     xlab = "Estimated tsd", ylab = "Expected tsd",
     main = "Estimated vs Expected tsd",
     xlim = c(range(c(tsds, tsds_e))),
     ylim = c(range(c(tsds, tsds_e))))

mu <- results_list$mean[grep("mu_tmu", rownames(results_list))]
mu_e <- results_list$expected[grep("mu_tmu", rownames(results_list))]

boxplot(mu,
       main = "Estimated mu_tmu",
       ylab = "Estimated mu_tmu")
abline(h = mu_e[1], col = "blue", lwd = 2, lty = 3)
