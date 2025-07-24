########### Thermal limits of marine species based on occurrence data ##########
# July of 2025
# Authors: Silas Principe, Richard McElreath, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Tests for model 2 ##############################

library(rethinking)
source("functions/sim_dataset.R")
set.seed(2025)

# Settings
stan_model_version <- 2
n_species <- 5 # More than 5 species is failing in this version
n_pts <- 100
mu_priors <- c(15, 20, 25)
sd_priors <- c(2, 5)
priors_grid <- expand.grid(mu_priors = mu_priors, sd_priors = sd_priors)

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
    tomax = 1,                                   # Max occupancy prob per species
    site_min = -Inf,                             # Minimum of site (for truncation)
    site_max = Inf                               # Maximum of site (for truncation)
)

get_dataset_stats(dataset_1)

# Run test 1
m1_data <- prepare_data_stan(dataset_1)

par(mfrow = c(2,3))
for (i in seq_len(nrow(priors_grid))) {

    m1_data$P_mu <- priors_grid$mu_priors[i]
    m1_data$P_sd <- priors_grid$sd_priors[i]

    m1 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m1_data, rstan_out = FALSE)

    precis(m1, 2)
    m1_results <- extract_precis(m1, dataset_1)

    l <- with(m1_results[m1_results$what == "tmu",], lm(expected ~ mean))
    with(m1_results[m1_results$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected)),
          main = paste0("MU", m1_data$P_mu, "/ SD ", m1_data$P_sd, " - rsq: ", round(summary(l)$adj.r.squared, 2))))
    abline(v = 20, h = 20, lty = 2)
}


# Extreme conditions - truncated ------
dataset_2 <- sim_dataset(
    N = n_pts * n_species,                       # Number of surveys
    n_sp = n_species,                            # Number of species
    p = rbeta(1, 2, 2),                          # Detection probability
    mu_tmu = 20,                                 # Global mean for tmu
    sigma_tmu = 3,                               # Global sd for tmu
    mu_tsd = 5,                                  # Global mean for tsd
    sigma_tsd = 1,                               # Global sd for tsd 
    tomax = 1,                                   # Max occupancy prob per species
    site_min = -Inf,                             # Minimum of site (for truncation)
    site_max = 25                                # Maximum of site (for truncation)
)

get_dataset_stats(dataset_2)

# Run test 2
m2_data <- prepare_data_stan(dataset_2)

par(mfrow = c(2,3))
for (i in seq_len(nrow(priors_grid))) {

    m2_data$P_mu <- priors_grid$mu_priors[i]
    m2_data$P_sd <- priors_grid$sd_priors[i]

    m2 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m2_data, rstan_out = FALSE)

    precis(m2, 2)
    m2_results <- extract_precis(m2, dataset_1)

    l <- with(m2_results[m2_results$what == "tmu",], lm(expected ~ mean))
    with(m2_results[m2_results$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected)),
          main = paste0("MU", m2_data$P_mu, "/ SD ", m2_data$P_sd, " - rsq: ", round(summary(l)$adj.r.squared, 2))))
    abline(v = 20, h = 20, lty = 2)
}


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
    tomax = 1,                                   # Max occupancy prob per species
    site_min = -Inf,                             # Minimum of site (for truncation)
    site_max = Inf,                              # Maximum of site (for truncation)
    n_presence = 50,                             # Number of presences
    n_absence = 50                               # Number of absences
)

get_dataset_stats(dataset_3)

# Run test 3
n_abs <- seq(50, 0, -5)
grid_results <- vector("list", nrow(priors_grid))

for (i in seq_len(nrow(priors_grid))) {

    pg_main <- paste0("PrMU: ", priors_grid$mu_priors[i], " PrSD: ", priors_grid$sd_priors[i], "\n")
    results_test3 <- vector("list", length(n_abs))
    dataset_mod <- dataset_3
    for (k in seq_along(n_abs)) {
        groups <- split(dataset_mod$dataset, dataset_mod$dataset$sid)
        groups <- lapply(groups, \(x){
            p <- x[x$y == 1,]
            a <- x[x$y == 0,]
            if (nrow(a) > n_abs[k]) {
                a <- a[sample(nrow(a), n_abs[k]), ]
            }
            rbind(p, a)
        })
        dataset_mod$dataset <- as.data.frame(do.call("rbind", groups))
        m3_data <- prepare_data_stan(dataset_mod)
        m3_data$P_mu <- priors_grid$mu_priors[i]
        m3_data$P_sd <- priors_grid$sd_priors[i]
        m3 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m3_data, rstan_out = FALSE)

        results_test3[[k]] <- extract_precis(m3, dataset_3)
    }
    par(mfrow = c(3, 4))
    for (k in seq_along(results_test3)) {
        l <- with(results_test3[[k]][results_test3[[k]]$what == "tmu",], lm(expected ~ mean))
        with(results_test3[[k]][results_test3[[k]]$what == "tmu",],
        plot(expected, mean, pch = 19, col = "#0c76c2",
            xlab = "Expected", ylab = "Predicted",
            xlim = range(c(mean, expected)), ylim = range(c(mean, expected)),
            main = ifelse(k == 1,
                          paste0(pg_main, n_abs[k], " - rsq: ", round(summary(l)$adj.r.squared, 2)),
                          paste0(n_abs[k], " - rsq: ", round(summary(l)$adj.r.squared, 2)))))
        abline(v = 20, h = 20, lty = 2)
    }

    grid_results[[i]] <- results_test3
}


# Real species ------
fish_data <- read.csv("data/fish_pa_data.csv")
head(fish_data)
fish_data <- fish_data[,c("species", "sst", "presence")]
colnames(fish_data) <- c("sid", "sst", "y")
fish_data$species <- fish_data$sid
fish_data$sid <- as.integer(as.factor(fish_data$sid))
# This model is failing when many species, restricting to 5
n_presences <- aggregate(fish_data$y, list(fish_data$sid), sum)[,2]
n_absences <- aggregate(fish_data$y, list(fish_data$sid), length)[,2] - n_presences
fish_data <- fish_data[fish_data$sid %in% unique(fish_data$sid)[n_presences > 5],]
fish_data <- fish_data[fish_data$sid %in% sample(unique(fish_data$sid), 5),]

# Run test 4
m4_data <- prepare_data_stan(list(dataset = fish_data))

par(mfrow = c(2,3))
for (i in seq_len(nrow(priors_grid))) {

    m4_data$P_mu <- priors_grid$mu_priors[i]
    m4_data$P_sd <- priors_grid$sd_priors[i]

    m4 <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m4_data, rstan_out = FALSE)

    precis(m4, 2)
    m4_results <- extract_precis(m4, dataset_1)

    l <- with(m4_results[m4_results$what == "tmu",], lm(expected ~ mean))
    with(m4_results[m4_results$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected)),
          main = paste0("MU", m4_data$P_mu, "/ SD ", m4_data$P_sd, " - rsq: ", round(summary(l)$adj.r.squared, 2))))
    abline(v = 20, h = 20, lty = 2)
}
