########### Thermal limits of marine species based on occurrence data ##########
# July of 2025
# Authors: Silas Principe, Richard McElreath, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Tests for model 7 ##############################

library(rethinking)
set.seed(2025)
source("functions/sim_dataset.R")

# Settings
stan_model_version <- 7
n_species <- 30
n_pts <- 100

N <- n_pts * n_species                       # Number of surveys
n_sp <- n_species                            # Number of species
p <- rbeta(1, 2, 2)                          # Detection probability
mu_tmu <- 20                                 # Global mean for tmu
sigma_tmu <- 3                               # Global sd for tmu
mu_tsd <- 5                                  # Global mean for tsd
sigma_tsd <- 1                               # Global sd for tsd 
tomax <- rbeta(n_species, 5, 1)              # Max occupancy prob per species - you can also pass a fixed value
site_min <- -Inf                             # Minimum of site (for truncation)
site_max <- Inf                              # Maximum of site (for truncation)
n_presence <- NULL                           # Number of presences
n_absence <- NULL                            # Number of absences
perc_sst_inc <- 0                            # Percentage of increase in SST SD

tmu <- round(rnorm(n_sp, mu_tmu, sigma_tmu), 1)
tsd <- round(abs(rnorm(n_sp, mu_tsd, sigma_tsd)), 1)

# Species IDs
sid <- rep(seq_len(n_sp), each = N / n_sp)

# Simulate SST for each survey
if (site_min != -Inf || site_max != Inf) {
    sst <- truncnorm::rtruncnorm(n = N, a = site_min, b = site_max, mean = tmu[sid], sd = tsd[sid])
} else {
    sst <- rnorm(N, mean = tmu[sid], sd = (tsd[sid] + (tsd[sid] * perc_sst_inc)))
}

gen_dataset <- function(sst, tmu, tsd, sid, tomax, p) {
    # Calculate occupancy probability (Gaussian suitability * tomax)
    q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2) * tomax[sid]

    # Simulate detections
    y <- numeric(N)
    occupancy <- numeric(N)
    for (i in seq_len(N)) {
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

# Dataset 1, normal
dataset_1 <- gen_dataset(sst, tmu, tsd, sid, tomax, p)

# Dataset 2-7 : varying p, fixed tomax
p_vector <- seq(0.5, 1, by = 0.1)
for (i in seq_along(p_vector)) {
    eval(parse(text = paste0("dataset_", i+1, " <- gen_dataset(sst, tmu, tsd, sid, rep(1, max(sid)), p=", p_vector[i], ")")))
}

# Dataset 8-13 : fixed p, varying tomax
tomax_vector <- seq(0.5, 1, by = 0.1)
for (i in seq_along(p_vector)) {
    eval(parse(text = paste0("dataset_", i+7, " <- gen_dataset(sst, tmu, tsd, sid, tomax = rep(tomax_vector[", i,"], max(sid)), p=1)")))
}

# Test models
results <- vector(mode = "list", length = 13)
results_m5 <- results

for (i in seq_len(13)) {
    message("Model ", i, "\n\n")
    eval(parse(text = paste0("sel_dataset <- dataset_", i)))
    m_data <- prepare_data_stan(sel_dataset)
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m_data, rstan_out = FALSE)
    m_5 <- cstan(file = paste0("codes/model5.stan"), data = m_data, rstan_out = FALSE)

    results[[i]] <- extract_precis(m, sel_dataset, transform = T, change_names = T)
    results_m5[[i]] <- extract_precis(m_5, sel_dataset, transform = T, change_names = T)
}

par(mfrow = c(4, 4))
titles <- c("All random", paste("Fixed tomax | p=", p_vector), paste("Fixed p | tomax=", tomax_vector))
for (i in seq_len(13)) {
    with(results[[i]][results[[i]]$what == "tmu",],
     plot(expected, mean, pch = 19, col = "#0c76c2",
          xlab = "Expected", ylab = "Predicted",
          xlim = range(c(mean, expected)), ylim = range(c(mean, expected)),
          main = titles[i]))
    abline(v = 20, h = 20, lty = 2)
}

ba_plot <- function(expected, predicted, main = NULL) {
    mean_me <- (expected + predicted) / 2
    diff_me <- expected - predicted
    loa <- mean(diff_me) + c(-1.96, 1.96) * sd(diff_me)

    plot(mean_me, diff_me,
        pch = 19, cex = 0.8,
        xlab = "Mean",
        ylab = "Difference expected - predicted",
        main = main, ylim = range(c(loa, diff_me)))

    abline(h = mean(diff_me), col = "#a11552", lwd = 2)
    abline(h = loa, col = "#1559ac", lty = 2)

    return(invisible(NULL))
}
par(mfrow = c(4, 4))
for (i in seq_len(13)) {
    with(results[[i]][results[[i]]$what == "tmu",],
     ba_plot(expected, mean, main = titles[i]))
}

lapply(results, \(dt) {
    dt <- dt[dt$what == "tmu",]
    round(summary(lm(expected ~ mean, data = dt))$adj.r.squared, 2)
})

abs_delta <- unlist(lapply(2:7, \(x) abs(results[[x]]$delta[results[[x]]$what == "tmu"])), use.names = F)
abs_delta_5 <- unlist(lapply(2:7, \(x) abs(results_m5[[x]]$delta[results_m5[[x]]$what == "tmu"])), use.names = F)
par(mfrow = c(1, 1))
plot(abs_delta, pch = 19, col = "#005e99", ylim = c(0, max(c(abs_delta, abs_delta_5)) + 0.5))
points(abs_delta_5, col = "#005e99")
for (i in 1:6) abline(v = seq(30, 180, by = 30)[i], lty = 2)
for (i in 1:6) text(x = seq(30, 180, by = 30)[i] - 15, y = max(c(abs_delta, abs_delta_5)) + 0.5, p_vector[i])
legend("right",
       legend = c("Model 7", "Model 5"),
       col = c("#005e99", "#005e99"),
       pch = c(19, 1))


# Dataset 14-18 : fixed p, varying tomax
n_total <- seq(100, 10, by = -20)
for (i in seq_along(n_total)) {
    sampled_ds <- dataset_1
    sampled_ds$dataset <- sampled_ds$dataset |>
        dplyr::group_by(sid) |>
        dplyr::slice_sample(n = n_total[i])
    eval(parse(text = paste0(
        "dataset_", i+13, "<- sampled_ds"
    )))
}

for (i in 14:18) {
    message("Model ", i, "\n\n")
    eval(parse(text = paste0("sel_dataset <- dataset_", i)))
    m_data <- prepare_data_stan(sel_dataset)
    m <- cstan(file = paste0("codes/model", stan_model_version, ".stan"), data = m_data, rstan_out = FALSE)
    m_5 <- cstan(file = paste0("codes/model5.stan"), data = m_data, rstan_out = FALSE)

    results[[i]] <- extract_precis(m, sel_dataset, transform = T, change_names = T)
    results_m5[[i]] <- extract_precis(m_5, sel_dataset, transform = T, change_names = T)
}

abs_delta <- unlist(lapply(14:18, \(x) abs(results[[x]]$delta[results[[x]]$what == "tmu"])), use.names = F)
abs_delta_5 <- unlist(lapply(14:18, \(x) abs(results_m5[[x]]$delta[results_m5[[x]]$what == "tmu"])), use.names = F)
mean_delta <- unlist(lapply(14:18, \(x) mean(abs(results[[x]]$delta[results[[x]]$what == "tmu"]))), use.names = F)
mean_delta_5 <- unlist(lapply(14:18, \(x) mean(abs(results_m5[[x]]$delta[results_m5[[x]]$what == "tmu"]))), use.names = F)
par(mfrow = c(1, 1))
plot(abs_delta, pch = 19, col = "#005e99", ylim = c(0, max(c(abs_delta, abs_delta_5)) + 0.5))
points(abs_delta_5, col = "#005e99")
for (i in 1:6) abline(v = seq(30, 180, by = 30)[i], lty = 2)
for (i in 1:6) text(x = seq(30, 180, by = 30)[i] - 15, y = max(c(abs_delta, abs_delta_5)) + 0.5, n_total[i])
legend("right",
       legend = c("Model 7", "Model 5"),
       col = c("#005e99", "#005e99"),
       pch = c(19, 1))
for (i in 1:6) lines(y = c(mean_delta[i], mean_delta[i]), x = c(seq(0, 150, by = 30)[i], seq(30, 180, by = 30)[i]))
for (i in 1:6) lines(y = c(mean_delta_5[i], mean_delta_5[i]), x = c(seq(0, 150, by = 30)[i], seq(30, 180, by = 30)[i]), lty = 2)
