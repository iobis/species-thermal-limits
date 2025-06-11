# Load packages
set.seed(2023)
library(rethinking)
library(dplyr)
library(ggplot2)
# For truncated normal
if (!require(truncnorm, quietly = T)) {
    install.packages("truncnorm")
}

# Data from https://doi.org/10.1002/ece3.6407
thermal_limits <- bind_rows(
    read.csv("data-raw/ece36407-sup-0002-datas2.csv"),
    read.csv("data-raw/ece36407-sup-0003-datas3.csv")
)

# Reef Life Survey
rls <- bind_rows(
    lapply(1:3, \(x) read.csv(paste0("data-raw/rls_globalreeffish_", x, ".csv")))
)

# Check species in common
sp_across <- thermal_limits$Species[thermal_limits$Species %in% rls$Taxon]

# Filter
thermal_limits <- thermal_limits |>
    filter(Species %in% sp_across)

get_pa <- function(species) {

    base_r <- terra::rast(res = 0.1)

    present <- rls |>
        filter(Taxon == species) |>
        select(x = SiteLong, y = SiteLat)

    absent <- rls |>
        filter(Taxon != species) |>
        select(x = SiteLong, y = SiteLat)
    
    base_r[terra::cellFromXY(base_r, absent)] <- 0
    base_r[terra::cellFromXY(base_r, present)] <- 1

    od <- as.data.frame(base_r, xy = T)
    colnames(od) <- c("decimalLongitude", "decimalLatitude", "presence")
    od$species <- species
    return(od)
}

pa_datasets <- lapply(unique(thermal_limits$Species), get_pa)

# Extract temperature
sst <- terra::rast("data-raw/thetao_baseline_depthsurf_mean.tif")

pa_datasets <- lapply(pa_datasets, \(ds) {
    temp <- terra::extract(sst, ds[,1:2], ID = F)
    ds$sst <- temp[,1]
    ds <- ds[!is.na(ds$sst),]
    ds
})

to_remove <- unlist(lapply(seq_along(pa_datasets), \(id) {
    n_pres <- sum(pa_datasets[[id]]$presence)
    if (n_pres < 10) {
        TRUE
    } else {
        FALSE
    }
}))

pa_datasets[to_remove] <- NULL

final_ds <- bind_rows(pa_datasets)

final_ds |> group_by(species, presence) |> count() |> filter(presence == 1)

# Run models
dat <- list(
    N = nrow(final_ds),
    N_spp = length(unique(final_ds$species)),
    sid = as.integer(as.factor(final_ds$species)),
    sst = final_ds$sst,
    y = final_ds$presence
)
m <- cstan(file = paste0("codes/model5.stan"), data = dat, rstan_out = FALSE)

prec_res <- precis(m, 2)

# Plot species
plot_sp <- function(sid) {

    sp_name <- levels(as.factor(final_ds$species))[sid]

    mean_val <- prec_res$mean[row.names(prec_res) == paste0("tmu[", sid, "]")]
    sd_val <- prec_res$mean[row.names(prec_res) == paste0("tsd[", sid, "]")]
    
    x <- seq(mean_val - 4 * sd_val, mean_val + 4 * sd_val, length.out = 100)
    y <- dnorm(x, mean = mean_val, sd = sd_val)

    plot(x, y, type = "l", lwd = 2, col = "blue",
        main = sp_name, xlab = "Value", ylab = "Density")

    sp_ex_dat <- final_ds[final_ds$species == sp_name,]
    data_mean <- mean(sp_ex_dat$sst[sp_ex_dat$presence == 1])

    abline(v = data_mean, col = "orange")

    sp_dat <- thermal_limits[thermal_limits$Species == sp_name,]

    abline(v = sp_dat$mean_CT_max, col = "#ff0073")
    abline(v = sp_dat$bo_sst_mean, col = "#343434")
    abline(v = sp_dat$bo_sst_q95, col = "#00a1d2")
    legend("topright",                          # Position
       legend = c("Data mean", "CT max", "Paper mean", "Paper95"),  # Labels
       col = c("orange", "#ff0073", "#343434", "#00a1d2"),               # Colors
       lty = 1, 
       lwd = 2, 
       bty = "n") 

    cli::cli_h1(sp_name)
    print(sp_dat[,c("mean_CT_max", "bo_sst_q5", "bo_sst_mean", "bo_sst_q95")])
    return(invisible(NULL))
}

# Use the ID to plot a specific species
plot_sp(1)
plot_sp(2)
plot_sp(3)
plot_sp(4)
prec_res1 <- prec_res

m6 <- cstan(file = paste0("codes/model6.stan"), data = dat, rstan_out = FALSE)

prec_res <- precis(m6, 2)
plot_sp(7)

# With less absences
species_list <- split(final_ds, final_ds$species)

# Function to process each species group
subsample_function <- function(species_df) {
  present <- species_df[species_df$presence == 1, ]
  absent <- species_df[species_df$presence == 0, ]
  
  # Subsample up to 50 zeros
  if (nrow(absent) > 10) {
    absent <- absent[sample(nrow(absent), 10), ]
  }
  
  # Combine and return
  rbind(present, absent)
}

# Apply the function to each species
result_list <- lapply(species_list, subsample_function)
df_subsampled <- do.call(rbind, result_list)

# Run models
dat <- list(
    N = nrow(df_subsampled),
    N_spp = length(unique(df_subsampled$species)),
    sid = as.integer(as.factor(df_subsampled$species)),
    sst = df_subsampled$sst,
    y = df_subsampled$presence
)
m2 <- cstan(file = paste0("codes/model5.stan"), data = dat, rstan_out = FALSE)

prec_res2 <- precis(m2, 2)

plot_sp2 <- function(sid) {

    sp_name <- levels(as.factor(final_ds$species))[sid]

    mean_val <- prec_res2$mean[row.names(prec_res2) == paste0("tmu[", sid, "]")]
    sd_val <- prec_res2$mean[row.names(prec_res2) == paste0("tsd[", sid, "]")]
    
    x <- seq(mean_val - 4 * sd_val, mean_val + 4 * sd_val, length.out = 100)
    y <- dnorm(x, mean = mean_val, sd = sd_val)

    plot(x, y, type = "l", lwd = 2, col = "blue",
        main = sp_name, xlab = "Value", ylab = "Density")

    sp_ex_dat <- final_ds[final_ds$species == sp_name,]
    data_mean <- mean(sp_ex_dat$sst[sp_ex_dat$presence == 1])

    abline(v = data_mean, col = "orange")

    sp_dat <- thermal_limits[thermal_limits$Species == sp_name,]

    abline(v = sp_dat$mean_CT_max, col = "#ff0073")
    abline(v = sp_dat$bo_sst_mean, col = "#343434")
    abline(v = sp_dat$bo_sst_q95, col = "#00a1d2")
    legend("topright",                          # Position
       legend = c("Data mean", "CT max", "Paper mean", "Paper95"),  # Labels
       col = c("orange", "#ff0073", "#343434", "#00a1d2"),               # Colors
       lty = 1, 
       lwd = 2, 
       bty = "n") 

    cli::cli_h1(sp_name)
    print(sp_dat[,c("mean_CT_max", "bo_sst_q5", "bo_sst_mean", "bo_sst_q95")])
    return(invisible(NULL))
}

# Use the ID to plot a specific species
plot_sp2(1)
plot_sp(1)

plot_sp2(2)
plot_sp(2)

plot_sp2(3)
plot_sp(3)
