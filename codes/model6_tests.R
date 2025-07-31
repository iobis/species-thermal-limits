########### Thermal limits of marine species based on occurrence data ##########
# July of 2025
# Authors: Silas Principe, Richard McElreath, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Tests for model 6 ##############################

library(rethinking)
source("functions/sim_dataset.R")
set.seed(2025)

# Settings
stan_model_version <- 6
n_species <- 30
n_pts <- 100

# Normal conditions ------
# We pass a fixed value for tomax, since the model 1 does not account for that
# By passing N=1000 and n_sp = 10 we are creating species with 100 records
dataset_1 <- sim_dataset(
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
m1_data <- prepare_data_stan(dataset_1)
m1 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m1_data, rstan_out = FALSE)

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

# Compare with model 5
m1_vs_5 <- cstan(file = "codes/model5.stan", data = m1_data, rstan_out = FALSE)

precis(m1_vs_5, 2)
m1_vs_5_results <- extract_precis(m1_vs_5, dataset_1)

# Mean
tmu_m1 <- m1_results$mean[m1_results$what == "tmu"]
tmu_m1_vs_5 <- m1_vs_5_results$mean[m1_vs_5_results$what == "tmu"]
d1_n_presence <- aggregate(dataset_1$dataset$y, list(dataset_1$dataset$sid), sum)[,2]

plot(tmu_m1[order(d1_n_presence)], pch = 19, col = "#005e99", ylim = range(c(tmu_m1, tmu_m1_vs_5)))
points(tmu_m1_vs_5[order(d1_n_presence)], col = "#005e99")
abline(h = dataset_1$mu_tmu, lty = 2)
legend("topright",
       legend = c("Model 6", "Model 5", "Global mu"),
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
       legend = c("Model 6", "Model 5"),
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
       legend = c("Model 6", "Model 5", "Global mu"),
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
       legend = c("Model 6", "Model 5"),
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
