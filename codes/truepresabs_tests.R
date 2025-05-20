################# Species Thermal Limits - modeling experiments ################
##################### Tests with true presence/absence data ####################
# Dec of 2024
# Author: Silas Principe, Pieter Provoost, Richard McElreath
# Contact: s.principe@unesco.org

# Dataset was prepared using the code on codes/data-processing/get_fish_pa.R
# It is based on data from the Reef Life Survey and the GASPAR project

set.seed(2023)
library(arrow)
library(dplyr)
library(rethinking)

# Create a folder to hold the results
res_f <- "results/simulations"
fs::dir_create(res_f)

fish_data <- arrow::read_parquet("data/fish_pa_data.parquet")
fish_data <- fish_data %>%
    filter(!is.na(surface_temperature))

occ_n <- table(fish_data$scientificname, fish_data$presence)
sp_good <- occ_n[occ_n[,2] >= 30,]

fish_data <- fish_data %>%
    filter(scientificname %in% row.names(sp_good))

# Remove absence data that is too far from presence
remove_absence <- function(presence, h3) {
    pres_h3 <- h3[presence == 1]
    # H3 at resolution 7 has a hex area of ~ 5km2, thus ~ 2.7km diagonal
    # 25 cells is thus ~ 55km appart
    near_disks <- h3jsr::get_disk(pres_h3, ring_size = 50)
    near_disks <- unique(unlist(near_disks))
    h3 %in% near_disks
}

fish_data  <- fish_data %>%
    group_by(scientificname) %>%
    mutate(to_keep = remove_absence(presence, h3_7)) %>%
    filter(to_keep) %>%
    select(-to_keep)

table(fish_data$scientificname, fish_data$presence)

# Sample 20 species for first test
sel_sp <- sample(unique(fish_data$scientificname), 20)
fish_data <- fish_data %>%
    filter(scientificname %in% sel_sp)
table(fish_data$scientificname, fish_data$presence)

#### TEMPORARY
to_remove <- function(presence, target = 500) {
    which_0 <- which(presence == 0)
    which_0 <- sample(which_0, target)
    to_r <- rep(FALSE, length(presence))
    to_r[presence == 0] <- TRUE
    to_r[which_0] <- FALSE
    to_r
}
fish_data <- fish_data %>%
    group_by(scientificname) %>%
    mutate(remove = to_remove(presence)) %>%
    filter(!remove) %>%
    select(-remove)

####

fish_summs <- fish_data %>%
    filter(presence == 1) %>%
    summarise(
        mean_s = mean(surface_temperature),
        sd_s = sd(surface_temperature),
        mean_b = mean(bottom_temperature),
        sd_b = sd(bottom_temperature)
    )

# Create object for modelling
dat <- list(
    N = nrow(fish_data),
    N_spp = length(unique(fish_data$scientificname)),
    sid = as.integer(as.factor(fish_data$scientificname)),
    sst = fish_data$surface_temperature,
    y = fish_data$presence
)

str(dat)

m_true1 <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

prec_res <- precis(m_true1, 2)
prec_res$expected <- c(fish_summs$mean_s, fish_summs$sd_s, 0)

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

write.csv(prec_res, file.path(res_f, "t2_1_true_pa_data.csv"), row.names = F)

# Investigate higher deviances
means_v <- prec_res[grepl("tmu", prec_res$what),]
means_v$delta <- means_v$expected - means_v$mean
means_v$species <- levels(as.factor(fish_data$scientificname))
means_v <- means_v[order(abs(means_v$delta), decreasing = T),]
means_v

worst_species <- means_v$species[1:5]
best_species <- means_v$species[15:20]

fish_data %>%
    filter(scientificname %in% worst_species) %>%
    group_by(scientificname, presence) %>%
    summarise(
        mean_sst = mean(surface_temperature),
        sd_sst = sd(surface_temperature)
    )

wrld <- rnaturalearth::ne_countries(returnclass = "sf")

fish_pts <- fish_data %>%
    filter(scientificname %in% worst_species)
fish_xy <- h3jsr::cell_to_point(fish_pts$h3_7)
fish_xy <- sf::st_as_sf(fish_xy)
fish_xy$species <- fish_pts$scientificname
fish_xy$presence <- fish_pts$presence

library(ggplot2)
ggplot() +
    geom_sf(data = wrld, fill = "grey80", color = "grey80") +
    geom_sf(data = fish_xy, aes(color = as.factor(presence)), alpha = .8) +
    theme_classic() +
    facet_wrap(~ species)


fish_data %>%
    filter(scientificname %in% worst_species) %>%
    group_by(scientificname, presence) %>%
    count()

fish_data %>%
    filter(scientificname %in% worst_species) %>%
    group_by(scientificname, presence) %>%
    summarise(mean_sst = mean(surface_temperature), sd_sst = sd(surface_temperature))

fish_data %>%
    filter(scientificname %in% best_species) %>%
    group_by(scientificname, presence) %>%
    ggplot() +
        geom_boxplot(aes(x = as.factor(presence), y = surface_temperature)) +
        facet_wrap(~scientificname)

fish_data %>%
    filter(scientificname %in% worst_species) %>%
    group_by(scientificname, presence) %>%
    ggplot() +
        geom_boxplot(aes(x = as.factor(presence), y = surface_temperature)) +
        facet_wrap(~scientificname)



fake_data <- data.frame(
    presence = c(rep(1, 30), rep(0, 30)),
    sst = truncnorm::rtruncnorm(60, mean = 27, sd = 1, b = 31),
    species = "Species A"
)

fake_data <- rbind(fake_data,
            data.frame(
    presence = c(rep(1, 30), rep(0, 30)),
    sst = truncnorm::rtruncnorm(60, mean = 28, sd = 1, b = 31),
    species = "Species B"
))

par(mfrow = c(1,1))
boxplot(fake_data$sst ~fake_data$presence+fake_data$species)


# Create object for modelling
dat <- list(
    N = nrow(fake_data),
    N_spp = length(unique(fake_data$species)),
    sid = as.integer(as.factor(fake_data$specie)),
    sst = fake_data$sst,
    y = fake_data$presence
)

str(dat)

m_true2 <- cstan(file = "codes/model4.stan", data = dat, rstan_out = FALSE)

prec_res2 <- precis(m_true2, 2)
prec_res2$expected <- c(27, 28, 1, 1, 0)
prec_res2

# Simulate data with lower upper bound for the suitability to
#see if same data