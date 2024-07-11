# Generate H3 indexed temperature dataset

# Load packages ----
library(terra)
library(h3r)
library(dplyr)
library(sf)
library(sfarrow)
library(arrow)

# Load/prepare environmental data ----
thetao <- rast(c("data-raw/thetao_baseline_depthsurf_mean.tif",
                 "data-raw/thetao_baseline_depthmean_mean.tif"))

names(thetao) <- c("surface", "bottom")

thetao <- as.data.frame(thetao, xy = T)

# Convert to H3 ----
thetao_h3 <- latLngToCell(lat = thetao$y, lng = thetao$x, 7)

thetao_c <- cbind(thetao, h3_7 = thetao_h3)

thetao_a <- thetao_c %>%
  select(h3_7, surface, bottom) %>%
  group_by(h3_7) %>%
  summarise(across(1:2, ~ mean(.x, na.rm = TRUE)))

head(thetao_a)

h3_xy <- cellToLatLng(thetao_a$h3_7)

thetao_a <- cbind(thetao_a, x = h3_xy$lng, y = h3_xy$lat)

thetao_a$surface <- round(thetao_a$surface, 2)
thetao_a$bottom <- round(thetao_a$bottom, 2)
thetao_a <- st_as_sf(thetao_a, coords = c("x", "y"), crs = st_crs(4326))
head(thetao_a)

st_write_parquet(thetao_a, "data/temperature_grid.parquet")

rm(list = ls())