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
    estimated_values <- data.frame(
        q_025 = exp(apply(results_opt_2[[i]]$f, 2:3, quantile, .25))[,1],
        median = exp(apply(results_opt_2[[i]]$f, 2:3, median)[,1]),
        q_075 = exp(apply(results_opt_2[[i]]$f, 2:3, quantile, .75))[,1]
    )
    n_pres <- aggregate(datasets_opt_2[[i]]$data$y, list(datasets_opt_2[[i]]$data$sid), sum)[,2]
    plot(estimated_values$median ~ true_values, main = paste("Optimum =", range_opt[i]),
         pch = 19, col = "#0c5094", xlim = c(0, 40))
    arrows(x0 = true_values, y0 = estimated_values$q_025, x1 = true_values, y1 = estimated_values$q_075, 
       angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)
    abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values$median ~ true_values), col = "grey50")
    text(x = true_values+1, y = estimated_values$median+1, labels = n_pres, cex = 1.5, col = "grey80")
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


# Try with equal number of presences, a low max_cov (ranges close to optimum)
n_pts <- 30

acanthurus <- rfishbase::species()
acanthurus <- acanthurus[acanthurus$Genus == "Acanthurus",]

acnth_estimates <- rfishbase::estimate()
acnth_estimates <- acnth_estimates[acnth_estimates$SpecCode %in% acanthurus$SpecCode,]

acnth_pref <- na.omit(acnth_estimates$TempPrefMean)
acnth_pref_m <- median(acnth_pref)
acnth_pref_sd <- sd(acnth_pref)

opt_t <- acnth_pref_m
max_cov <- 0.1
n_species <- 10

dataset_4 <- sim_dataset_phylo(
        N = n_pts * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(opt_t),
        max_cov = max_cov,
        rho_val = 1,
        equal_n_pres = TRUE
    )
m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = dataset_4$data, rstan_out = FALSE)
results_4 <- extract.samples(m)

par(mfrow = c(1,1))
true_values <- exp(dataset_4$raw$f[, 1])
estimated_values <- data.frame(
    q_025 = exp(apply(results_4$f, 2:3, quantile, .25))[,1],
    median = exp(apply(results_4$f, 2:3, median)[,1]),
    q_075 = exp(apply(results_4$f, 2:3, quantile, .75))[,1]
)
n_pres <- aggregate(dataset_4$data$y, list(dataset_4$data$sid), sum)[,2]
plot(estimated_values$median ~ true_values, main = paste("Optimum =", opt_t, "- equal N"),
        pch = 19, col = "#0c5094", xlim = c(15, 31), ylim = c(15, 31))
arrows(x0 = true_values, y0 = estimated_values$q_025, x1 = true_values, y1 = estimated_values$q_075, 
    angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)
#abline(v = c(5, 30), lty = 2)
abline(lm(estimated_values$median ~ true_values), col = "grey50")
text(x = true_values+1, y = estimated_values$median+1, labels = n_pres, cex = 1.5, col = "grey80")





#### New tests 6 Nov
opt_t <- 25
n_species <- 10

# Test 1 - Increasing max covariance (the higher, the higher the spread)
max_cov <- seq(0.1, 2.1, by = 0.4)

max_cov_results <- max_cov_datasets <- vector("list", length(max_cov))

for (i in seq_along(max_cov)) {
    message("\nFitting model ", i, " out of ", length(max_cov), " ============\n\n")
    max_cov_datasets[[i]] <- sim_dataset_phylo(
        N = n_pts * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(opt_t),
        max_cov = max_cov[i],
        rho_val = 1,
        equal_n_pres = TRUE
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = max_cov_datasets[[i]]$data, rstan_out = FALSE)
    max_cov_results[[i]] <- extract.samples(m)
}

par(mfrow = c(2,3))
for (i in seq_along(max_cov_datasets)) {
    true_values <- exp(max_cov_datasets[[i]]$raw$f[, 1])
    estimated_values <- data.frame(
        q_025 = exp(apply(max_cov_results[[i]]$f, 2:3, quantile, .25))[,1],
        median = exp(apply(max_cov_results[[i]]$f, 2:3, median)[,1]),
        q_075 = exp(apply(max_cov_results[[i]]$f, 2:3, quantile, .75))[,1]
    )
    n_pres <- aggregate(max_cov_datasets[[i]]$data$y, list(max_cov_datasets[[i]]$data$sid), sum)[,2]
    plot(estimated_values$median ~ true_values, 
         main = paste("Maximum covariance =", max_cov[i], "- equal N - optimum =", opt_t),
         pch = 19, col = "#0c5094")#, xlim = c(10, 40), ylim = c(10, 40))
    arrows(x0 = true_values, y0 = estimated_values$q_025, x1 = true_values, y1 = estimated_values$q_075, 
        angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)
    #abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values$median ~ true_values), col = "grey50")
    text(x = true_values+1, y = estimated_values$median+1, labels = n_pres, cex = 1.5, col = "grey80")
}


# Test 2 - Increasing number of points
max_cov <- 0.1
opt_t <- 25
n_species <- 10

records <- seq(30, 280, by = 50)
records_results <- records_datasets <- vector("list", length(records))

for (i in seq_along(records)) {
    message("\nFitting model ", i, " out of ", length(records), " ============\n\n")
    records_datasets[[i]] <- sim_dataset_phylo(
        N = records[i] * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(opt_t),
        max_cov = max_cov,
        rho_val = 1,
        equal_n_pres = TRUE
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = records_datasets[[i]]$data, rstan_out = FALSE)
    records_results[[i]] <- extract.samples(m)
}

par(mfrow = c(2,3))
for (i in seq_along(records_datasets)) {
    true_values <- exp(records_datasets[[i]]$raw$f[, 1])
    estimated_values <- data.frame(
        q_025 = exp(apply(records_results[[i]]$f, 2:3, quantile, .25))[,1],
        median = exp(apply(records_results[[i]]$f, 2:3, median)[,1]),
        q_075 = exp(apply(records_results[[i]]$f, 2:3, quantile, .75))[,1]
    )
    n_pres <- aggregate(records_datasets[[i]]$data$y, list(records_datasets[[i]]$data$sid), sum)[,2]
    plot(estimated_values$median ~ true_values, 
         main = paste("N records =", records[i], "- equal N - optimum =", opt_t),
         pch = 19, col = "#0c5094")#, xlim = c(10, 40), ylim = c(10, 40))
    arrows(x0 = true_values, y0 = estimated_values$q_025, x1 = true_values, y1 = estimated_values$q_075, 
        angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)
    #abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values$median ~ true_values), col = "grey50")
    text(x = true_values+1, y = estimated_values$median+1, labels = n_pres, cex = 1.5, col = "grey80")
}









# Test 3 - Different amounts of presence/absence
max_cov <- 0.1
opt_t <- 25
n_species <- 10

# Create a dataset with a large number of records, which will be subsampled
pa_dataset <- sim_dataset_phylo(
    N = 200 * n_species,
    n_sp = n_species,
    p = 0.8, 
    global_optimum = log(opt_t),
    max_cov = max_cov,
    rho_val = 1,
    equal_n_pres = TRUE
)

ab_n <- seq(60, 0, -10)
pr_n <- 60

pa_sub_datasets <- lapply(ab_n, \(n) {
    edited_ds <- pa_dataset
    edited_ds$data |>
        dplyr::group_by()
})

pa_results <- vector("list", length(pa_sub_datasets))

for (i in seq_along(records)) {
    message("\nFitting model ", i, " out of ", length(records), " ============\n\n")
    records_datasets[[i]] <- sim_dataset_phylo(
        N = records[i] * n_species,                       # Number of surveys
        n_sp = n_species,                            # Number of species
        p = 0.8, 
        global_optimum = log(opt_t),
        max_cov = max_cov,
        rho_val = 1,
        equal_n_pres = TRUE
    )
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = records_datasets[[i]]$data, rstan_out = FALSE)
    records_results[[i]] <- extract.samples(m)
}

par(mfrow = c(2,3))
for (i in seq_along(records_datasets)) {
    true_values <- exp(records_datasets[[i]]$raw$f[, 1])
    estimated_values <- data.frame(
        q_025 = exp(apply(records_results[[i]]$f, 2:3, quantile, .25))[,1],
        median = exp(apply(records_results[[i]]$f, 2:3, median)[,1]),
        q_075 = exp(apply(records_results[[i]]$f, 2:3, quantile, .75))[,1]
    )
    n_pres <- aggregate(records_datasets[[i]]$data$y, list(records_datasets[[i]]$data$sid), sum)[,2]
    plot(estimated_values$median ~ true_values, 
         main = paste("N records =", records[i], "- equal N - optimum =", opt_t),
         pch = 19, col = "#0c5094")#, xlim = c(10, 40), ylim = c(10, 40))
    arrows(x0 = true_values, y0 = estimated_values$q_025, x1 = true_values, y1 = estimated_values$q_075, 
        angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)
    #abline(v = c(5, 30), lty = 2)
    abline(lm(estimated_values$median ~ true_values), col = "grey50")
    text(x = true_values+1, y = estimated_values$median+1, labels = n_pres, cex = 1.5, col = "grey80")
}