################# Species Thermal Limits - modeling experiments ################
# August of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Check Reef Life Survey and GASPAR dataset for presence-absence data

# Load packages ----
library(worrms)
library(dplyr)
library(tidyr)
library(arrow)
library(terra)

# Load data
gaspar <- read.csv("data-raw/gaspar_dataset.csv")
gaspar <- gaspar %>%
  select(decimalLongitude = Longitude, decimalLatitude = Latitude,
         species = Species, abundance = Abundance)

rls <- do.call("rbind", 
               lapply(list.files("data-raw/", pattern = "rls_", full.names = T),
                      read.csv))
rls <- rls %>%
  select(decimalLongitude = SiteLong, decimalLatitude = SiteLat,
         species = Taxon, abundance = Total)


fish_abund <- bind_rows(gaspar, rls)

eg <- worrms::wm_records_names("Acanthurus chirurgus")
eg <- eg[[1]]
eg[,] <- NA
fish_tax <- lapply(unique(fish_abund$species), function(x) {
  fr <- try(worrms::wm_records_names(x))
  if (inherits(fr, "try-error")) {
    return(eg)
  } else {
    fr <- fr[[1]]
    if (any(grepl("^accepted", fr$status))) {
      fr <- fr[fr$status == "accepted",]
    } else {
      rec <- wm_record(fr$valid_AphiaID[1])
      if (any(grepl("^accepted", rec$status))) {
        fr <- rec[rec$status == "accepted",]
      } else {
        fr <- rec[1,]
      }
    }
    return(fr)
  }
})

fish_tax <- do.call("rbind", fish_tax)
fish_tax$species <- unique(fish_abund$species)

fish_abund <- left_join(fish_abund, fish_tax)
fish_abund <- fish_abund[fish_abund$rank == "Species",]

# Add H3
fish_abund <- fish_abund %>% 
  filter(decimalLongitude <= 180 & decimalLongitude >= -180) %>%
  filter(decimalLatitude <= 90 & decimalLatitude >= -90)
fish_abund$h3_7 <- h3jsr::point_to_cell(fish_abund[,c("decimalLongitude", "decimalLatitude")], res = 7)

fish_abund_per_area <- fish_abund %>%
  group_by(scientificname, AphiaID, class, family, h3_7) %>%
  summarise(abundance = sum(abundance))

# Create presence absence data.frame
fish_pa <- fish_abund_per_area %>%
  pivot_wider(names_from = h3_7, values_from = abundance)

fish_pa[is.na(fish_pa)] <- 0

fish_pa <- fish_pa %>% 
  pivot_longer(5:ncol(.), names_to = "h3_7", values_to = "presence")

fish_pa$presence[fish_pa$presence > 0] <- 1

# Add environmental information
# Download SST data from Bio-ORACLE
download.file("https://erddap.bio-oracle.org/erddap/griddap/thetao_baseline_2000_2019_depthsurf.nc?thetao_mean%5B(2000-01-01):1:(2010-01-01T00:00:00Z)%5D%5B(-89.975):1:(89.975)%5D%5B(-179.975):1:(179.975)%5D",
"data-raw/thetao_baseline_depthsurf_mean.nc")
download.file("https://erddap.bio-oracle.org/erddap/griddap/thetao_baseline_2000_2019_depthmean.nc?thetao_mean%5B(2000-01-01):1:(2010-01-01T00:00:00Z)%5D%5B(-89.975):1:(89.975)%5D%5B(-179.975):1:(179.975)%5D",
"data-raw/thetao_baseline_depthmean_mean.nc")

sst <- rast("data-raw/thetao_baseline_depthsurf_mean.nc")
smdt <- rast("data-raw/thetao_baseline_depthmean_mean.nc")

sst <- mean(sst)
smdt <- mean(smdt)

writeRaster(sst, "data-raw/thetao_baseline_depthsurf_mean.tif", overwrite = T)
writeRaster(smdt, "data-raw/thetao_baseline_depthmean_mean.tif", overwrite = T)

file.remove(c("data-raw/thetao_baseline_depthsurf_mean.nc", "data-raw/thetao_baseline_depthmean_mean.nc"))

# Extract
data_points <- h3jsr::cell_to_point(unique(fish_pa$h3_7), simple = F)

# Function to fill NAs by neighbor
get_nearest <- function(x, y, rlayer) {
    cell <- cellFromXY(rlayer, cbind(x, y))
    look <- matrix(c(rep(1, 12), 0, rep(1, 12)), nrow = 5, ncol = 5, byrow = T)
    adj <- adjacent(rlayer, cell, directions = look)
    vals <- rlayer[as.vector(adj)]
    vals <- vals[,1]
    if (all(is.na(vals))) {
        return(NA)
    } else if (length(na.omit(vals)) == 1) {
        return(na.omit(vals))
    } else {
        adj_ok <- adj[which(!is.na(vals))]
        closer <- nearest(
          vect(cbind(x, y), crs = "EPSG:4326"),
          vect(
            xyFromCell(rlayer, adj_ok),
            crs = "EPSG:4326")
          )
        vals[!is.na(vals)][closer$to_id]
    }
}

sst_data <- terra::extract(sst, data_points, ID = F)[,1]

sst_data[is.na(sst_data)] <- unlist(lapply(which(is.na(sst_data)), function(id){
  pt <- data_points[id,]
  pt <- sf::st_coordinates(pt)
  get_nearest(pt[,1], pt[,2], sst)
}))

smdt_data <- terra::extract(smdt, data_points, ID = F)[,1]

smdt_data[is.na(smdt_data)] <- unlist(lapply(which(is.na(smdt_data)), function(id){
  pt <- data_points[id,]
  pt <- sf::st_coordinates(pt)
  get_nearest(pt[,1], pt[,2], smdt)
}))

sst_full <- data.frame(
  h3_7 = data_points$h3_address,
  surface_temperature = sst_data,
  bottom_temperature = smdt_data
)

fish_pa <- left_join(fish_pa, sst_full)
fish_abund_per_area <- left_join(fish_abund_per_area, sst_full)

# Save as a parquet file
arrow::write_parquet(fish_pa, "data/fish_pa_data.parquet")
arrow::write_parquet(fish_abund_per_area, "data/fish_abund_data.parquet")
