########### Thermal limits of marine species based on occurrence data ##########
# July of 2025
# Authors: Silas Principe, Richard McElreath, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Tests for model 7 ##############################

library(rethinking)
source("functions/sim_dataset.R")
set.seed(2025)

# Settings
stan_model_version <- 8
n_species <- c(10, 30, 100, 200)
n_pts <- 30

results <- vector("list", length = length(n_species))
datasets <- results
times <- rep(NA, length(datasets))

for (i in seq_along(n_species)) {

    message("Testing with ", n_species[i], " species.\n\n")

    datasets[[i]] <- sim_dataset_phylo(
        N = n_pts * n_species[i],                       # Number of surveys
        n_sp = n_species[i],                            # Number of species
        p = 0.8, global_optimum = log(15), max_covariance = 0.5, rho = 0.2
    )

    tictoc::tic()
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = datasets[[i]]$data, rstan_out = FALSE)
    times[i] <- tictoc::toc(log = T, quiet = T)$callback_msg

    results[[i]] <- extract.samples(m)

}

par(mfrow = c(1,3))
for (i in 1:3) {
    f_est <- apply(results[[i]]$f, 2:3, mean)
    plot(datasets[[i]]$raw$f[, 1], f_est[, 1], xlab = "simulated", ylab = "estimated")
    abline(a = 0, b = 1, lty = 2)
}
for (i in 1:3) {
    f_est <- apply(results[[i]]$f, 2:3, mean)
    plot(datasets[[i]]$raw$f[, 2], f_est[, 2], xlab = "simulated", ylab = "estimated")
    abline(a = 0, b = 1, lty = 2)
}
for (i in 1:3) {
    f_est <- apply(results[[i]]$f, 2:3, mean)
    plot(datasets[[i]]$raw$f[, 3], f_est[, 3], xlab = "simulated", ylab = "estimated")
    abline(a = 0, b = 1, lty = 2)
}


# Test effect of different optimums
# Site min is 5
# Site max is 30
range_opt <- seq(5, 25, by = 5)
results_opt <- vector("list", length = length(range_opt))
datasets_opt <- results_opt
n_species <- 10
n_pts <- 30

for (i in seq_along(range_opt)) {
    message("\nFitting model ", i, " out of ", length(range_opt), "============\n\n")
    datasets_opt[[i]] <- sim_dataset_phylo(
        N = n_pts * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(range_opt[i])
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = datasets_opt[[i]]$data, rstan_out = FALSE)
    results_opt[[i]] <- extract.samples(m)
}

par(mfrow = c(2,3))
for (i in seq_along(datasets_opt)) {
    true_values <- exp(datasets_opt[[i]]$raw$f[, 1])
    estimated_values <- exp(apply(results_opt[[i]]$f, 2:3, mean)[,1])
    n_pres <- aggregate(datasets_opt[[i]]$data$y, list(datasets_opt[[i]]$data$sid), sum)[,2]
    plot(estimated_values ~ true_values, main = paste("Range opt =", range_opt[i]),
         pch = 19, col = "#0c5094", xlim = c(0, 40))
    abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values ~ true_values), col = "grey50")
    text(x = true_values, y = estimated_values+1, labels = n_pres, cex = 1.5)
}

# Try with equal number of presences
range_opt <- seq(5, 25, by = 5)
results_opt_2 <- vector("list", length = length(range_opt))
datasets_opt_2 <- results_opt_2

for (i in seq_along(range_opt)) {
    message("\nFitting model ", i, " out of ", length(range_opt), " ============\n\n")
    datasets_opt_2[[i]] <- sim_dataset_phylo(
        N = n_pts * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(range_opt[i]),
        equal_n_pres = TRUE
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = datasets_opt_2[[i]]$data, rstan_out = FALSE)
    results_opt_2[[i]] <- extract.samples(m)
}

par(mfrow = c(2,3))
for (i in seq_along(datasets_opt_2)) {
    true_values <- exp(datasets_opt_2[[i]]$raw$f[, 1])
    estimated_values <- exp(apply(results_opt_2[[i]]$f, 2:3, mean)[,1])
    n_pres <- aggregate(datasets_opt_2[[i]]$data$y, list(datasets_opt_2[[i]]$data$sid), sum)[,2]
    plot(estimated_values ~ true_values, main = paste("Optimum =", range_opt[i]),
         pch = 19, col = "#0c5094", xlim = c(0, 40))
    abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values ~ true_values), col = "grey50")
    text(x = true_values, y = estimated_values+1, labels = n_pres, cex = 1.5)
}

# Try with equal number of presences and now with a fixed opt, but varying maximum covariance
range_max_cov <- seq(1, 3, by = 0.5)
results_cov_2 <- vector("list", length = length(range_max_cov))
datasets_cov_2 <- results_cov_2

for (i in seq_along(range_max_cov)) {
    message("\nFitting model ", i, " out of ", length(range_max_cov), " ============\n\n")
    datasets_cov_2[[i]] <- sim_dataset_phylo(
        N = n_pts * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(20),
        equal_n_pres = TRUE,
        max_covariance = range_max_cov[i]
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = datasets_cov_2[[i]]$data, rstan_out = FALSE)
    results_cov_2[[i]] <- extract.samples(m)
}

par(mfrow = c(2,3))
for (i in seq_along(datasets_cov_2)) {
    true_values <- exp(datasets_cov_2[[i]]$raw$f[, 1])
    estimated_values <- exp(apply(results_cov_2[[i]]$f, 2:3, mean)[,1])
    n_pres <- aggregate(datasets_cov_2[[i]]$data$y, list(datasets_cov_2[[i]]$data$sid), sum)[,2]
    plot(estimated_values ~ true_values, main = paste("Maximum cov =", range_max_cov[i]),
         pch = 19, col = "#0c5094", xlim = c(0, 40))
    abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values ~ true_values), col = "grey50")
    text(x = true_values, y = estimated_values+1, labels = n_pres, cex = 1.5)
}

# Try with equal number of presences and now with a fixed opt, fixed maximum covariance, but varying rho
range_rho <- seq(0.5, 2, by = 0.5)
results_rho <- vector("list", length = length(range_rho))
datasets_rho <- results_rho

for (i in seq_along(range_max_cov)) {
    message("\nFitting model ", i, " out of ", length(range_max_cov), " ============\n\n")
    datasets_rho[[i]] <- sim_dataset_phylo(
        N = n_pts * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(20),
        equal_n_pres = TRUE,
        max_covariance = 3,
        rho_val = range_rho[i]
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = datasets_rho[[i]]$data, rstan_out = FALSE)
    results_rho[[i]] <- extract.samples(m)
}

par(mfrow = c(2,2))
for (i in seq_along(datasets_rho)) {
    true_values <- exp(datasets_rho[[i]]$raw$f[, 1])
    estimated_values <- exp(apply(results_rho[[i]]$f, 2:3, mean)[,1])
    n_pres <- aggregate(datasets_rho[[i]]$data$y, list(datasets_rho[[i]]$data$sid), sum)[,2]
    plot(estimated_values ~ true_values, main = paste("rho value =", range_rho[i]),
         pch = 19, col = "#0c5094", xlim = c(0, 40))
    abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values ~ true_values), col = "grey50")
    text(x = true_values, y = estimated_values+1, labels = n_pres, cex = 1.5)
}















# Normal conditions ------
# We pass a fixed value for tomax, since the model 1 does not account for that
# By passing N=1000 and n_sp = 10 we are creating species with 100 records
dataset_1 <- sim_dataset_phylo(
    N = n_pts * n_species,                       # Number of surveys
    n_sp = n_species,                            # Number of species
    p = rbeta(1, 2, 2),                          # Detection probability
    mu_tmu = 20,                                 # Global mean for tmu
    sigma_tmu = 3,                               # Global sd for tmu
    mu_tsd = 5,                                  # Global mean for tsd
    sigma_tsd = 1,                               # Global sd for tsd 
    tomax = rbeta(n_species, 5, 1),              # Max occupancy prob per species (variable in this version)
    site_min = -Inf,                             # Minimum of site (for truncation)
    site_max = Inf                               # Maximum of site (for truncation)
)

get_dataset_stats(dataset_1)

# Run test 1
m1 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = dataset_1, rstan_out = FALSE)

precis(m1, 2)
m1_results <- extract_precis(m1, dataset_1, transform = T, change_names = T)



















with(m1_results[m1_results$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected))))
abline(v = 20, h = 20, lty = 2)
summary(with(m1_results[m1_results$what == "tmu",], lm(expected ~ mean)))

with(m1_results[m1_results$what == "tomax",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted", main = "tomax",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected))))
summary(with(m1_results[m1_results$what == "tomax",], lm(expected ~ mean)))

with(m1_results[m1_results$what == "tmu",],
     plot(y = delta, x = dataset_1$tmu, xlab = "tmu", ylab = "Delta", pch = 19))
abline(h = 0, lty = 2)
abline(v = 20, lty = 2, col = "blue")

# Compare with model 5
m1_vs_5 <- cstan(file = "codes/model5.stan", data = m1_data, rstan_out = FALSE)

precis(m1_vs_5, 2)
m1_vs_5_results <- extract_precis(m1_vs_5, dataset_1, transform = T, change_names = T)

# Mean
tmu_m1 <- m1_results$mean[m1_results$what == "tmu"]
tmu_m1_vs_5 <- m1_vs_5_results$mean[m1_vs_5_results$what == "tmu"]
d1_n_presence <- aggregate(dataset_1$dataset$y, list(dataset_1$dataset$sid), sum)[,2]

plot(tmu_m1[order(d1_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m1, tmu_m1_vs_5)))
points(tmu_m1_vs_5[order(d1_n_presence)], col = "#005e99")
abline(h = dataset_1$mu_tmu, lty = 2)
legend("topright",
       legend = c("Model 7", "Model 5", "Global mu"),
       col = c("#005e99", "#005e99", "black"),
       pch = c(19, 1, NA),
       lty = c(NA, NA, 2))

# Delta
tmu_m1 <- m1_results$delta[m1_results$what == "tmu"]
tmu_m1_vs_5 <- m1_vs_5_results$delta[m1_vs_5_results$what == "tmu"]

plot(tmu_m1[order(d1_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m1, tmu_m1_vs_5)))
points(tmu_m1_vs_5[order(d1_n_presence)], col = "#005e99")
abline(h = 0, lty = 2)
legend("topright",
       legend = c("Model 7", "Model 5"),
       col = c("#005e99", "#005e99"),
       pch = c(19, 1))


# Extreme conditions - truncated ------
dataset_2 <- sim_dataset(
    N = n_pts * n_species,                       # Number of surveys
    n_sp = n_species,                            # Number of species
    p = rbeta(1, 2, 2),                          # Detection probability
    mu_tmu = 20,                                 # Global mean for tmu
    sigma_tmu = 3,                               # Global sd for tmu
    mu_tsd = 5,                                  # Global mean for tsd
    sigma_tsd = 1,                               # Global sd for tsd 
    tomax = rbeta(n_species, 5, 1),              # Max occupancy prob per species (variable in this version)
    site_min = -Inf,                             # Minimum of site (for truncation)
    site_max = 26                                # Maximum of site (for truncation)
)

get_dataset_stats(dataset_2)

# Run test 2
m2_data <- prepare_data_stan(dataset_2)
m2 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m2_data, rstan_out = FALSE)

precis(m2, 2)
m2_results <- extract_precis(m2, dataset_2, transform = T, change_names = T)

with(m2_results[m2_results$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected))))
abline(v = 20, h = 20, lty = 2)
summary(with(m2_results[m2_results$what == "tmu",], lm(expected ~ mean)))

# Compare with model 5
m2_vs_5 <- cstan(file = "codes/model5.stan", data = m2_data, rstan_out = FALSE)

precis(m2_vs_5, 2)
m2_vs_5_results <- extract_precis(m2_vs_5, dataset_2, transform = T, change_names = T)

# Mean
tmu_m2 <- m2_results$mean[m2_results$what == "tmu"]
tmu_m2_vs_5 <- m2_vs_5_results$mean[m2_vs_5_results$what == "tmu"]
d2_n_presence <- aggregate(dataset_2$dataset$y, list(dataset_2$dataset$sid), sum)[,2]

plot(tmu_m2[order(d2_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m2, tmu_m2_vs_5)))
points(tmu_m2_vs_5[order(d2_n_presence)], col = "#005e99")
abline(h = dataset_2$mu_tmu, lty = 2)
legend("topright",
       legend = c("Model 7", "Model 5", "Global mu"),
       col = c("#005e99", "#005e99", "black"),
       pch = c(19, 1, NA),
       lty = c(NA, NA, 2))

# Delta
tmu_m2 <- m2_results$delta[m2_results$what == "tmu"]
tmu_m2_vs_5 <- m2_vs_5_results$delta[m2_vs_5_results$what == "tmu"]

plot(tmu_m2[order(d2_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m2, tmu_m2_vs_5)))
points(tmu_m2_vs_5[order(d2_n_presence)], col = "#005e99")
abline(h = 0, lty = 2)
legend("topright",
       legend = c("Model 7", "Model 5"),
       col = c("#005e99", "#005e99"),
       pch = c(19, 1))


# Varying number of absences ------
# We pass a fixed value for tomax, since the model 1 does not account for that
# By passing N=1000 and n_sp = 10 we are creating species with 100 records
dataset_3 <- sim_dataset(
    N = n_pts * n_species,                       # Number of surveys
    n_sp = n_species,                            # Number of species
    p = rbeta(1, 2, 2),                          # Detection probability
    mu_tmu = 20,                                 # Global mean for tmu
    sigma_tmu = 3,                               # Global sd for tmu
    mu_tsd = 5,                                  # Global mean for tsd
    sigma_tsd = 1,                               # Global sd for tsd 
    tomax = rbeta(n_species, 5, 1),              # Max occupancy prob per species (variable in this version)
    site_min = -Inf,                             # Minimum of site (for truncation)
    site_max = Inf,                              # Maximum of site (for truncation)
    n_presence = 50,                             # Number of presences
    n_absence = 50                               # Number of absences
)

get_dataset_stats(dataset_3)

# Run test 3
n_abs <- seq(50, 0, -5)
results_test3 <- vector("list", length(n_abs))
stan_data <- vector("list", length(n_abs))
dataset_mod <- dataset_3
for (i in seq_along(n_abs)) {
    message("\nRunning ", i, " out of ", length(n_abs), "\n")
    groups <- split(dataset_mod$dataset, dataset_mod$dataset$sid)
    groups <- lapply(groups, \(x){
        p <- x[x$y == 1,]
        a <- x[x$y == 0,]
        if (nrow(a) > n_abs[i]) {
            a <- a[sample(nrow(a), n_abs[i]), ]
        }
        rbind(p, a)
    })
    dataset_mod$dataset <- as.data.frame(do.call("rbind", groups))
    m3_data <- prepare_data_stan(dataset_mod)
    stan_data[[i]] <- m3_data
    m3 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m3_data, rstan_out = FALSE)

    results_test3[[i]] <- extract_precis(m3, dataset_3, transform = T, change_names = T)
}

par(mfrow = c(3, 4))
for (i in seq_along(results_test3)) {
    l <- with(results_test3[[i]][results_test3[[i]]$what == "tmu",], lm(expected ~ mean))
    with(results_test3[[i]][results_test3[[i]]$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected)),
          main = paste0(n_abs[i], " - rsq: ", round(summary(l)$adj.r.squared, 2))))
    abline(v = 20, h = 20, lty = 2)
}

# Compare with model 5
results_test3_m5 <- vector("list", length(stan_data))
for (i in seq_along(stan_data)) {
    m_temp <- cstan(file = "codes/model5.stan", data = stan_data[[i]], rstan_out = FALSE)
    results_test3_m5[[i]] <- extract_precis(m_temp, dataset_3, transform = T, change_names = T)
}

# Mean
par(mfrow = c(3, 4))
for (i in seq_along(stan_data)) {
    m3_vs_5_results <- results_test3_m5[[i]]
    m3_results <- results_test3[[i]]

    tmu_m3 <- m3_results$mean[m3_results$what == "tmu"]
    tmu_m3_vs_5 <- m3_vs_5_results$mean[m3_vs_5_results$what == "tmu"]
    d3_n_presence <- aggregate(dataset_3$dataset$y, list(dataset_3$dataset$sid), sum)[,2]

    plot(tmu_m3[order(d3_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m3, tmu_m3_vs_5)))
    points(tmu_m3_vs_5[order(d3_n_presence)], col = "#005e99")
    abline(h = dataset_3$mu_tmu, lty = 2)
    legend("topright",
        legend = c("Model 7", "Model 5", "Global mu"),
        col = c("#005e99", "#005e99", "black"),
        pch = c(19, 1, NA),
        lty = c(NA, NA, 2))
}

# Delta
par(mfrow = c(3, 4))
for (i in seq_along(stan_data)) {
    m3_vs_5_results <- results_test3_m5[[i]]
    m3_results <- results_test3[[i]]

    tmu_m3 <- m3_results$delta[m3_results$what == "tmu"]
    tmu_m3_vs_5 <- m3_vs_5_results$delta[m3_vs_5_results$what == "tmu"]
    d3_n_presence <- aggregate(dataset_3$dataset$y, list(dataset_3$dataset$sid), sum)[,2]

    plot(tmu_m3[order(d3_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m3, tmu_m3_vs_5)))
    points(tmu_m3_vs_5[order(d3_n_presence)], col = "#005e99")
    abline(h = dataset_3$mu_tmu, lty = 2)
    legend("topright",
        legend = c("Model 7", "Model 5", "Global mu"),
        col = c("#005e99", "#005e99", "black"),
        pch = c(19, 1, NA),
        lty = c(NA, NA, 2))
}


# Real species ------
fish_data <- read.csv("data/fish_pa_data.csv")
head(fish_data)
fish_data <- fish_data[,c("species", "sst", "presence")]
colnames(fish_data) <- c("sid", "sst", "y")
fish_data$species <- fish_data$sid
fish_data$sid <- as.integer(as.factor(fish_data$sid))

# Run test 4
m4_data <- prepare_data_stan(list(dataset = fish_data))
m4 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m4_data, rstan_out = FALSE)

precis(m4, 2)
m4_results <- extract_precis(m4, list(
    dataset = fish_data,
    tmu = aggregate(fish_data$sst[fish_data$y == 1], list(fish_data$sid[fish_data$y == 1]), mean)[,2],
    tsd = aggregate(fish_data$sst[fish_data$y == 1], list(fish_data$sid[fish_data$y == 1]), sd)[,2],
    p = 0, mu_tmu = NA, mu_tsd = NA, sigma_tmu = NA, sigma_tsd = NA
), transform = T, change_names = T)

par(mfrow = c(1,1))
with(m4_results[m4_results$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected))))
abline(v = 20, h = 20, lty = 2)
summary(with(m4_results[m4_results$what == "tmu",], lm(expected ~ mean)))




install.packages("fishtree")
library(fishtree)
library(ape)
phy <- fishtree_phylogeny(rank = "Acanthuridae")
phy
par(mfrow=c(2, 1))
plot(phy, show.tip.label = FALSE)
ltt.plot(phy)
