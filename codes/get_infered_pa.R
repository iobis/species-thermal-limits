# Get an infered presence-absence dataset

library(robis)
library(h3jsr)
library(dplyr)
library(terra)
coords <- c("decimalLongitude", "decimalLatitude")

target_species <- "Labrus viridis"

occ_data <- occurrence(target_species, absence = "include")
occ_data <- occ_data %>%
  select(aphiaID, species, family, familyid, decimalLongitude, decimalLatitude, absence) %>%
  rowwise() %>%
  mutate(h3_06 = h3jsr::point_to_cell(c(decimalLongitude, decimalLatitude), res = 6))

occ_area <- vect(occ_data, geom = coords, crs = "EPSG:4326")
occ_area <- terra::aggregate(terra::buffer(occ_area, width = 100000))

plot(occ_area)
points(occ_data[,coords], pch = 20, cex = .5)

occ_area_ext <- sf::st_as_text(sf::st_as_sfc(sf::st_bbox(occ_area)))

occ_family <- occurrence(taxonid = occ_data$familyid[1], absence = "include",
                         geometry = occ_area_ext)

occ_family <- occ_family %>%
  select(aphiaID, species, family, familyid, decimalLongitude, decimalLatitude, absence, date_mid) %>%
  mutate(date = as.POSIXct(date_mid / 1000, origin = "1970-01-01")) %>%
  mutate(year = lubridate::year(date), month = lubridate::month(date)) %>%
  filter(!is.na(date))

# Check where it was more sampled
occ_family$h3_06 <- h3jsr::point_to_cell(occ_family[,coords], res = 6)

occ_area_h3 <- h3jsr::polygon_to_cells(sf::st_as_sf(occ_area), res = 6)

occ_family_sel <- occ_family[occ_family$h3_06 %in% unlist(occ_area_h3),]

absences_family <- occ_family_sel %>%
  group_by(h3_06) %>%
  summarise(absence_f = all(absence)) %>%
  filter(absence_f)
# 0 here

occ_family_sample <- occ_family_sel %>%
  filter(aphiaID != occ_data$aphiaID[1]) %>%
  group_by(h3_06) %>%
  distinct(year, month) %>%
  summarise(effort = n())

non_valid <- h3jsr::get_ring(occ_data$h3_06, ring_size = 2)

occ_family_sample <- occ_family_sample %>%
  filter(!h3_06 %in% c(unlist(non_valid), occ_data$h3_06)) %>%
  filter(effort > 1)

absences <- sample(occ_family_sample$h3_06, 20, prob = occ_family_sample$effort)

absences_coords <- h3jsr::cell_to_point(absences)

lines(occ_area)
points(occ_data[,coords], pch = 20, cex = .5, col = "blue")
points(sf::st_coordinates(absences_coords), pch = 20, cex = .5, col = "red")

sst <- rast("../../mpa_europe/mpaeu_sdm/data/env/current/thetao_baseline_depthsurf_max.tif")
plot(sst)
