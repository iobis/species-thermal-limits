library(dplyr)

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
    if (n_pres < 1) {
        TRUE
    } else {
        FALSE
    }
}))

pa_datasets[to_remove] <- NULL

final_ds <- bind_rows(pa_datasets)

final_ds |> group_by(species, presence) |> count() |> filter(presence == 1)

write.csv(final_ds, "data/fish_pa_data.csv", row.names = F)
