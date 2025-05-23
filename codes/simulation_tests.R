################# Species Thermal Limits - modeling experiments ################
########################## Tests with simulated data ###########################
# Nov of 2024
# Author: Silas Principe, Pieter Provoost, Richard McElreath
# Contact: s.principe@unesco.org

set.seed(2023)
library(terra)
library(dplyr)
library(rethinking)
source("codes/simulate_species.R")

# Create a folder to hold the results and data
sim_f <- "data/simulated_data"
res_f <- "results/simulations"
fs::dir_create(sim_f)
fs::dir_create(res_f)


# TEST 1: basic occupancy model with thermal tolerance function with fixed prior

# Set environmental parameters (3 scenarios)
mean_site <- c(15, 20, 25)
sd_site <- c(5, 5, 5) # fixed SD for simplicity
cat(paste(paste("Site", seq_along(mean_site), ">>", (mean_site - 2*sd_site), "to", (mean_site + 2*sd_site)), collapse = "\n"))

# Set the number of species
n_species <- 40

# Set number of cells
n_cells <- 4000

# Set the means for each species
sel_x_species <- sample(15:28, n_species, replace = T)

# And also the variance
sel_xhat_species <- sample(1:4, n_species, replace = T)


# Simulate several species with different number
simulated_datasets <- lapply(1:(length(mean_site)*n_species), function(x) NULL)
counter <- 0
for (site in seq_along(mean_site)) {
    for (sp in seq_len(n_species)) {
        counter <- counter + 1
        sim_r <- sim_species(
            x_species = sel_x_species[sp],
            xhat_species = sel_xhat_species[sp],
            x_site = mean_site[site],
            xhat_site = sd_site[site],
            ncells = n_cells,
            # The upper and lower limits of environmental data were established
            # based on data from Bio-ORACLE v3.0. We considered the average
            # SST, although a test with the maximum could also be done
            site_min = -2,
            site_max = 31
        )
        sim_r <- as.data.frame(sim_r)
        sim_r$site <- paste("site", site)
        sim_r$species <- paste("species", sp)
        simulated_datasets[[counter]] <- sim_r
    }
}

simulated_datasets <- do.call("rbind", simulated_datasets)
simulated_datasets$species_site <- paste(simulated_datasets$species, simulated_datasets$site)

write.csv(simulated_datasets, file.path(sim_f, "simulation1.csv"), row.names = F)


# Check that none have no presences
# This should not be a problem when the simulation considers equal class sample
# as done in this case
summaries <- simulated_datasets |>
    dplyr::group_by(species, site) |>
    dplyr::summarise(
        presence = sum(sampled_occurrence),
        absence = length(sampled_occurrence) - sum(sampled_occurrence)
    )
if (any(summaries$presence == 0)) warning("No presence data!")
if (any(summaries$presence < 2)) warning("Sites with low number of presences")

# TEST 1.1 Equal number of presences and absences
# Make data object

dat <- list(
    N = nrow(simulated_datasets),
    N_spp = length(unique(simulated_datasets$species_site)),
    sid = as.integer(as.factor(simulated_datasets$species_site)),
    sst = simulated_datasets$surface,
    y = simulated_datasets$sampled_occurrence
)

spp_mu <- rep(NA, dat$N_spp)
spp_sd <- spp_mu
for (i in seq_len(dat$N_spp)) {
    suit <- simulated_datasets$surface_suitability[dat$sid==i]
    sst <- simulated_datasets$surface[dat$sid==i]
    spp_mu[i] <- sst[suit == max(suit)]
    sst_b <- simulated_datasets$surface[dat$sid==i & dat$y == 1]
    spp_sd[i] <- sd(sst_b)
}

m_sim1 <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

prec_res <- precis(m_sim1, 2)
prec_res$expected <- c(spp_mu, spp_sd, 0)

# Plot to check
par(mfrow = c(1,2))
plot(y = prec_res$mean[grepl("tmu", row.names(prec_res))],
     x = prec_res$expected[grepl("tmu", row.names(prec_res))],
     xlab = "Expected", ylab = "Predicted", main = "Mean", pch = 20, col = "#2c648e77")
abline(lm(prec_res$mean[grepl("tmu", row.names(prec_res))] ~ prec_res$expected[grepl("tmu", row.names(prec_res))]))
plot(y = prec_res$mean[grepl("tsd", row.names(prec_res))],
     x = prec_res$expected[grepl("tsd", row.names(prec_res))],
     xlab = "Expected", ylab = "Predicted", main = "SD", pch = 20, col = "#2c648e77")
abline(lm(prec_res$mean[grepl("tsd", row.names(prec_res))] ~ prec_res$expected[grepl("tsd", row.names(prec_res))]))

# Compare fits using a Bland-Altman plot
ba_plot_stats <- blandr::blandr.statistics(
    prec_res$expected[grepl("tmu", row.names(prec_res))],
    prec_res$mean[grepl("tmu", row.names(prec_res))]
)
blandr::blandr.plot.ggplot(ba_plot_stats) + 
    ggplot2::geom_vline(xintercept = 20) +
    ggplot2::geom_label(x = 21, y = -3.5, label = "Prior for mean")

# Save results
prec_res <- as.data.frame(prec_res)
prec_res$what <- row.names(prec_res)
row.names(prec_res) <- NULL
class(prec_res) <- "data.frame"

write.csv(prec_res, file.path(res_f, "t1_1_equal_number.csv"), row.names = F)


# TEST 1.2 Varying number of absences
target_n_absence <- rev(seq(0, 18, by = 2))
results_n_absence <- lapply(target_n_absence, function(x) NULL)

for (ab in seq_along(target_n_absence)) {
    # Make data object
    to_sample <- target_n_absence[ab]

    presences <- simulated_datasets %>%
        filter(sampled_occurrence == 1)

    absences <- simulated_datasets %>%
        filter(sampled_occurrence == 0) %>%
        group_by(species_site) %>%
        slice_sample(n = to_sample)

    reduced_dataset <- bind_rows(presences, absences)

    dat <- list(
        N = nrow(reduced_dataset),
        N_spp = length(unique(reduced_dataset$species_site)),
        sid = as.integer(as.factor(reduced_dataset$species_site)),
        sst = reduced_dataset$surface,
        y = reduced_dataset$sampled_occurrence
    )

    m_sim2 <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

    prec_res <- precis(m_sim2, 2)
    prec_res$expected <- c(spp_mu, spp_sd, 0)
    prec_res <- as.data.frame(prec_res)
    prec_res$n_absence <- to_sample

    prec_res$what <- row.names(prec_res)
    row.names(prec_res) <- NULL

    results_n_absence[[ab]] <- prec_res
}

results_n_absence <- do.call("rbind", results_n_absence)
results_n_absence <- as.data.frame(results_n_absence)

par(mfrow = c(2, ceiling(length(target_n_absence)/2)))
mean_results <- results_n_absence[grepl("tmu", results_n_absence$what),]
for (tp in target_n_absence) {
    fd <- mean_results[mean_results$n_absence == tp, ]
    mval <- ceiling(range(c(
        mean_results$mean, mean_results$expected
    )))
    plot(mean ~ expected, data = fd, main = paste(tp, "absences"), col = "#0066c077",
        xlim = mval, pch = 20)
    points(prec_res$mean[grepl("tmu", prec_res$what)] ~ prec_res$expected[grepl("tmu", prec_res$what)],
           col = "#e27f0e55", pch = 20)
}

write.csv(results_n_absence, file.path("results/simulations", "t1_2_varying_absences.csv"))


# TEST 1.3: Simulated dataset based on true surface

# Download SST data
# download.file("https://erddap.bio-oracle.org/erddap/griddap/thetao_baseline_2000_2019_depthsurf.nc?thetao_mean%5B(2000-01-01):1:(2010-01-01T00:00:00Z)%5D%5B(-89.975):1:(89.975)%5D%5B(-179.975):1:(179.975)%5D",
# "data-raw/thetao_baseline_depthsurf_mean.nc")
# sst <- rast("data-raw/thetao_baseline_depthsurf_mean.nc")
# sst <- mean(sst)
# writeRaster(sst, "data-raw/thetao_baseline_depthsurf_mean.tif", overwrite = T)
# file.remove("data-raw/thetao_baseline_depthsurf_mean.nc")

sst <- rast("data-raw/thetao_baseline_depthsurf_mean.tif")

# Crop to an area
cam <- ext(-86.308594, -58.095703, 12.46876, 32.916485)
sst_c <- crop(sst, cam)
plot(sst_c)
sst_c_vals <- values(sst_c)[,1]
sst_c_vals <- sst_c_vals[!is.na(sst_c_vals)]
length(sst_c_vals)
summary(sst_c_vals)

# Set the number of species
n_species <- 6

# Set number of cells
n_cells <- 4000

# Set the means for each species
sel_x_species <- seq(from = 23, to = 28)

# And also the variance
sel_xhat_species <- rep(1, n_species)


# Simulate several species with different number
simulated_ds_env <- lapply(1:n_species, function(x) NULL)
for (sp in seq_len(n_species)) {
    sim_r <- sim_species_env(
        x_species = sel_x_species[sp],
        xhat_species = sel_xhat_species[sp],
        site_values = sst_c_vals,
        ncells = n_cells
    )
    sim_r <- as.data.frame(sim_r)
    sim_r$species <- paste("species", sp)
    simulated_ds_env[[sp]] <- sim_r
}

simulated_ds_env <- do.call("rbind", simulated_ds_env)
table(simulated_ds_env$sampled_occurrence, simulated_ds_env$species)

# Create another dataset with equal number of absences by subsampling
simulated_ds_env_eq <- simulated_ds_env %>%
    group_by(species, sampled_occurrence) %>%
    sample_n(size = 150)
table(simulated_ds_env_eq$sampled_occurrence, simulated_ds_env_eq$species)

simulated_ds_env$type <- "original"
simulated_ds_env_eq$type <- "reduced"

simulated_ds_env <- rbind(simulated_ds_env, simulated_ds_env_eq)

write.csv(simulated_ds_env, file.path(sim_f, "simulation2.csv"), row.names = F)

# We will work only with the reduced for now
simulated_ds_env <- simulated_ds_env[simulated_ds_env$type == "reduced",]
for (i in unique(simulated_ds_env$species)) {
    spd <- simulated_ds_env[simulated_ds_env$species == i,]
    cat("Avg pres = ", mean(spd$surface[spd$true_occurrence == 1]),
    "Avg abs = ", mean(spd$surface[spd$true_occurrence == 0]), "\n")
    plot(density(spd$surface[spd$true_occurrence == 1]))
    lines(density(spd$surface[spd$true_occurrence == 0]), col = "red")
}

# Make data object
dat <- list(
    N = nrow(simulated_ds_env),
    N_spp = length(unique(simulated_ds_env$species)),
    sid = as.integer(as.factor(simulated_ds_env$species)),
    sst = simulated_ds_env$surface,
    y = simulated_ds_env$sampled_occurrence
)

# spp_mu <- rep(NA, dat$N_spp)
# spp_sd <- spp_mu
# for (i in seq_len(dat$N_spp)) {
#     suit <- simulated_ds_env$surface_suitability[dat$sid==i]
#     sst <- simulated_ds_env$surface[dat$sid==i]
#     spp_mu[i] <- sst[suit == max(suit)]
#     sst_b <- simulated_ds_env$surface[dat$sid==i & dat$y == 1]
#     spp_sd[i] <- sd(sst_b)
# }
spp_mu <- sel_x_species
spp_sd <- sel_xhat_species

m_sim3 <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

prec_res <- precis(m_sim3, 2)
prec_res$expected <- c(spp_mu, spp_sd, 0)

# Plot to check
par(mfrow = c(1,1))
plot(y = prec_res$mean[grepl("tmu", row.names(prec_res))],
     x = prec_res$expected[grepl("tmu", row.names(prec_res))],
     xlab = "Expected", ylab = "Predicted", main = "Mean", pch = 20, col = "#2c648e77")
abline(lm(prec_res$mean[grepl("tmu", row.names(prec_res))] ~ prec_res$expected[grepl("tmu", row.names(prec_res))]))

# Compare fits using a Bland-Altman plot
ba_plot_stats <- blandr::blandr.statistics(
    prec_res$expected[grepl("tmu", row.names(prec_res))],
    prec_res$mean[grepl("tmu", row.names(prec_res))]
)
blandr::blandr.plot.ggplot(ba_plot_stats) + 
    ggplot2::geom_vline(xintercept = 20) +
    ggplot2::geom_label(x = 21, y = -3.5, label = "Prior for mean")

# Save results
prec_res <- as.data.frame(prec_res)
prec_res$what <- row.names(prec_res)
row.names(prec_res) <- NULL
class(prec_res) <- "data.frame"

write.csv(prec_res, file.path(res_f, "t1_1_equal_number.csv"), row.names = F)

# TEST 1.4: Simulated dataset based on true surface with varying number of absences

# TEST 1.5: Simulated dataset based on true surface with varying bias