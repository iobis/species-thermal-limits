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

# Add H3
fish_abund$h3_7 <- h3jsr::point_to_cell(fish_abund[,c("decimalLongitude", "decimalLatitude")], res = 7)

fish_abund_pa <- fish_abund %>%
  group_by(scientificname, h3_7) %>%
  summarise(abundance = sum(abundance))


