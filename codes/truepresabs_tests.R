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

fish_summs <- fish_data %>%
    filter(presence == 1) %>%
    summarise(
        mean_s = mean(surface_temperature),
        sd_s = sd(surface_temperature),
        mean_b = mean(bottom_temperature),
        sd_b = sd(bottom_temperature)
    )

# Remove absence data that is too far from presence
remove_absence <- function(presence, h3) {
    pres_h3 <- h3[presence == 1]
    near_disks <- h3jsr::get_disk(pres_h3)
    near_disks <- unique(unlist(near_disks))
    h3 %in% near_disks
}

fish_data  <- fish_data %>%
    group_by(scientificname) %>%
    mutate(to_keep = remove_absence(presence, h3_7)) %>%
    filter(to_keep) %>%
    select(-to_keep)

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
