################# Species Thermal Limits - modeling experiments ################
########################## Tests with simulated data ###########################
# Nov of 2024
# Author: Silas Principe, Pieter Provoost, Richard McElreath
# Contact: s.principe@unesco.org

set.seed(2023)
library(dplyr)
library(rethinking)
source("codes/simulate_species.R")

# Create a folder to hold the results and data
fs::dir_create("data/simulated_data")
fs::dir_create("results/simulations")

# TEST 1: basic occupancy model with thermal tolerance function with fixed prior

# Set environmental parameters (5 scenarios)
mean_site <- c(15, 20)
sd_site <- c(10, 5) # fixed SD for simplicity
cat(paste(paste("Site", seq_along(mean_site), ">>", (mean_site - 2*sd_site), "to", (mean_site + 2*sd_site)), collapse = "\n"))

# Set the number of species
n_species <- 20

# Set number of cells
n_cells <- 1000

# Set the means for each species
sel_x_species <- sample(15:28, n_species, replace = T)#c(29, 27, 26, 20, 15)

# And also the variance
sel_xhat_species <- sample(1:4, n_species, replace = T)#c(2, 3, 1, 2, 4)


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

# Check that none have no presences
summaries <- simulated_datasets |>
    dplyr::group_by(species, site) |>
    dplyr::summarise(
        presence = sum(sampled_occurrence),
        absence = length(sampled_occurrence) - sum(sampled_occurrence)
    )
any(summaries$presence == 0)
any(summaries$presence < 2)

# TEST 1.1 Equal number of presences and absences
# Make data object
simulated_datasets$species_site <- paste(simulated_datasets$species, simulated_datasets$site)

# Only work with a subset. Still to be checked.
##simulated_datasets <- simulated_datasets[simulated_datasets$species_site %in% unique(simulated_datasets$species_site)[1:5],]

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
     xlab = "Expected", ylab = "Predicted", main = "Mean")
abline(lm(prec_res$mean[grepl("tmu", row.names(prec_res))] ~ prec_res$expected[grepl("tmu", row.names(prec_res))]))
plot(y = prec_res$mean[grepl("tsd", row.names(prec_res))],
     x = prec_res$expected[grepl("tsd", row.names(prec_res))],
     xlab = "Expected", ylab = "Predicted", main = "SD")
abline(lm(prec_res$mean[grepl("tsd", row.names(prec_res))] ~ prec_res$expected[grepl("tsd", row.names(prec_res))]))

prec_res <- as.data.frame(prec_res)
prec_res$what <- row.names(prec_res)
row.names(prec_res) <- NULL

write.csv(prec_res, file.path("results/simulations", "t1_1_equal_number.csv"))


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
    plot(mean ~ expected, data = fd, main = paste(tp, "absences"), col = "#0066c0",
        xlim = mval)
    points(prec_res$mean[grepl("tmu", prec_res$what)] ~ prec_res$expected[grepl("tmu", prec_res$what)],
           col = "#e27f0e", pch = 20)
}

write.csv(results_n_absence, file.path("results/simulations", "t1_2_varying_absences.csv"))



# TEST 1.3: Simulated dataset based on true surface

# TEST 1.4: Simulated dataset based on true surface with varying number of absences

# TEST 1.5: Simulated dataset based on true surface with varying error

# TEST 1.6: True dataset (fishes)

# Get datasets
rls <- do.call("rbind", lapply(list.files("data-raw", full.names = T, pattern = "rls"), read.csv))
rls$h3_6 <- h3jsr::point_to_cell(rls[,c("SiteLong", "SiteLat")], res = 6)

rls <- rls |>
    group_by(h3_6) |>
    distinct(Taxon)

pa_rls_matrix <- rls |>
    mutate(presence = 1) |>
    tidyr::pivot_wider(names_from = Taxon, values_from = presence)

pa_rls_matrix[is.na(pa_rls_matrix)] <- 0

pa_rls_ds <- pa_rls_matrix |>
    tidyr::pivot_longer(2:ncol(pa_rls_matrix), names_to = "species", values_to = "presence")

pa_rls_pts <- h3jsr::cell_to_point(pa_rls_ds$h3_6)

sst <- terra::rast("~/Research/mpa_europe/mpaeu_sdm/data/env/current/thetao_baseline_depthsurf_mean.tif")

sst_pts <- terra::extract(sst, terra::vect(pa_rls_pts), ID = F)[,1]

na_vals <- which(is.na(sst_pts))
na_pts <- terra::vect(pa_rls_pts)[na_vals,]

pb <- progress::progress_bar$new(total = length(na_vals))

for (i in seq_along(na_vals)) {
    pb$tick()
    tc <- terra::cellFromXY(sst, data.frame(
        x = terra::geom(na_pts[i,])[,c("x", "y")][1],
        y = terra::geom(na_pts[i,])[,c("x", "y")][2]
    ))
    ad <- terra::adjacent(sst, tc, directions = 8)
    vv <- terra::extract(sst, as.vector(ad))[,1]
    if (any(!is.na(vv))) {
        sst_pts[na_vals[i]] <- vv[!is.na(vv)][1]
    } 
}

pa_rls_ds$temperature <- sst_pts

pa_rls_ds <- pa_rls_ds[!is.na(pa_rls_ds$temperature),]

stats_sp <- pa_rls_ds |>
    group_by(species) |>
    summarise(presence_n = sum(presence), absence_n = sum(presence == 0))

stats_sp <- stats_sp |>
    ungroup() |>
    filter(presence_n > 20 & absence_n > 0)

pa_rls_ds <- pa_rls_ds[pa_rls_ds$species %in% stats_sp$species,]

dat <- list(
    N = nrow(pa_rls_ds),
    N_spp = length(unique(pa_rls_ds$species)),
    sid = as.integer(as.factor(pa_rls_ds$species)),
    sst = pa_rls_ds$temperature,
    y = pa_rls_ds$presence
)

spp_mu <- rep(NA, dat$N_spp)
spp_sd <- spp_mu
for (i in seq_len(dat$N_spp)) {
    spp_mu[i] <- mean(pa_rls_ds$temperature[dat$sid==i & dat$y == 1])
    spp_sd[i] <- sd(pa_rls_ds$temperature[dat$sid==i & dat$y == 1])
}

m_sim6 <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

prec_res <- precis(m_sim6, 2)
prec_res$expected <- c(spp_mu, spp_sd, 0)

# Plot to check
par(mfrow = c(1,2))
plot(y = prec_res$mean[grepl("tmu", row.names(prec_res))],
     x = prec_res$expected[grepl("tmu", row.names(prec_res))],
     xlab = "Expected", ylab = "Predicted", main = "Mean")
abline(lm(prec_res$mean[grepl("tmu", row.names(prec_res))] ~ prec_res$expected[grepl("tmu", row.names(prec_res))]))
plot(y = prec_res$mean[grepl("tsd", row.names(prec_res))],
     x = prec_res$expected[grepl("tsd", row.names(prec_res))],
     xlab = "Expected", ylab = "Predicted", main = "SD")
abline(lm(prec_res$mean[grepl("tsd", row.names(prec_res))] ~ prec_res$expected[grepl("tsd", row.names(prec_res))]))

prec_res <- as.data.frame(prec_res)
prec_res$what <- row.names(prec_res)
row.names(prec_res) <- NULL

write.csv(prec_res, file.path("results/simulations", "t1_6_true_species.csv"))






for (i in unique(simulated_datasets$species_site)) {
    d <- simulated_datasets[simulated_datasets$species_site == i,]
    plot(density(d$surface[d$sampled_occurrence == 1]))
    lines(density(d$surface[d$sampled_occurrence == 0]), col = "red")
}
