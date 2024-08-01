################# Species Thermal Limits - modeling experiments ################
####################### Meeting #2 - Exploring use cases #######################
# May of 2024
# Author: Silas Principe
# Contact: s.principe@unesco.org
#
# Check GlobTHERM available species

# To download GlobTHERM data: https://datadryad.org/stash/dataset/doi:10.5061/dryad.1cv08
# To download OBIS Full Export: https://obis.org/data/access/
# Alternatively, use OBIS R package `robis::occurrence()`

# Load packages ----
library(worrms)
library(dplyr)
library(tidyr)
library(arrow)
library(terra)

# Download files if needed
if (!file.exists("data-raw/GlobalTherm_upload_02_11_17.csv")) {
  res <- httr::GET("https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.1cv08/download",
                   httr::write_disk("data-raw/GlobTherm.zip", overwrite = TRUE))
  unzip("data-raw/GlobTherm.zip", exdir = "data-raw/")
  file.remove("data-raw/GlobTherm.zip")
}

if (!file.exists("data-raw/thetao_baseline_depthmean_mean.tif")) {
  ds <- c("thetao_baseline_2000_2019_depthsurf",
          "thetao_baseline_2000_2019_depthmean")
  for (dscode in ds) {
    # See: https://github.com/bio-oracle/biooracler
    bd <- biooracler::download_layers(dscode,
                                      variables = c("thetao_mean"),
                                      constraints = list(time = c('2001-01-01T00:00:00Z', '2010-01-01T00:00:00Z')))
    bd <- mean(bd)
    names(bd) <- "thetao_mean"
    
    writeRaster(bd, paste0("data-raw/", gsub("\\d+_\\d+_", "", dscode), "_mean.tif"))
    rm(bd)
  }
}

# Identify species from GlobTHERM that are marine ----
gb <- read.csv("data-raw/GlobalTherm_upload_02_11_17.csv")

gb <- gb %>%
  rowwise() %>%
  mutate(full_name = trimws(iconv(paste(Genus, Species), from = "ISO-8859-1", to = "UTF-8"),
                            whitespace = "[\\h\\v]"))

gb_species <- gb$full_name

gb_taxonomy <- obistools::match_taxa(gb_species, ask = F)

sp_list <- cbind(full_name = gb_species, gb_taxonomy)
sp_list <- sp_list[!is.na(sp_list$scientificName),]
sp_list <- sp_list[sp_list$match_type %in% c("exact", "phonetic"),]

sp_list_mar <- worrms::wm_record_(as.numeric(sp_list$acceptedNameUsageID))
sp_list_mar <- bind_rows(sp_list_mar)
sp_list_mar <- sp_list_mar[!is.na(sp_list_mar$isMarine),]
sp_list_mar <- sp_list_mar[sp_list_mar$isMarine == 1,]

sp_list <- sp_list[sp_list$scientificName %in% sp_list_mar$scientificname, ]

gb_selected <- gb %>%
  filter(full_name %in% sp_list$full_name) %>%
  select(full_name, Tmax, max_metric) %>%
  left_join(sp_list)

# Load OBIS/GBIF full export to get SST information
dbdata <- open_dataset("s3://obis-products/speciesgrids/h3_7/")

sst <- rast("data-raw/thetao_baseline_depthsurf_mean.tif")

dbdata_records <- dbdata %>%
  select(species, AphiaID, records, h3_07) %>%
  filter(AphiaID %in% gb_selected$acceptedNameUsageID) %>%
  collect()

dbdata_records_coord <- h3jsr::cell_to_point(dbdata_records$h3_07)
dbdata_records_coord <- sf::st_coordinates(dbdata_records_coord)
colnames(dbdata_records_coord) <- c("decimalLongitude", "decimalLatitude")

dbdata_records <- cbind(dbdata_records, dbdata_records_coord)

rec_sst <- extract(sst, dbdata_records[,c("decimalLongitude", "decimalLatitude")])

dbdata_records$sst <- rec_sst$thetao_mean

dbdata_records_sst <- dbdata_records %>%
  group_by(species, AphiaID) %>%
  filter(!is.na(sst)) %>%
  summarise(maxSST = max(sst, na.rm = T),
            meanSST = mean(sst, na.rm = T),
            sdSST = sd(sst, na.rm = T),
            medianSST = median(sst, na.rm = T),
            minSST = min(sst, na.rm = T),
            records = sum(records),
            n_cells = n())

gb_selected_f <- gb_selected %>%
  rename(AphiaID = acceptedNameUsageID) %>%
  mutate(AphiaID = as.numeric(AphiaID)) %>%
  left_join(dbdata_records_sst) %>%
  filter(!is.na(Tmax))

to_annotate <- gb_selected_f %>%
  filter(Tmax - maxSST > 20)

to_annotate_ft <- worrms::wm_record_(to_annotate$AphiaID)
to_annotate_ft <- bind_rows(to_annotate_ft)

to_annotate$phylum <- to_annotate_ft$phylum

ggplot(gb_selected_f) +
  geom_smooth(aes(x = Tmax, y = maxSST)) +
  geom_point(aes(x = Tmax, y = maxSST, color = sdSST), size = 2.5) +
  scale_color_viridis_c() +
  xlab("Maximum temperature (laboratory)") + ylab("Maximum temperature (occurrence data)") +
  geom_label(data = to_annotate, aes(x = Tmax + 6, y = maxSST + c(0.4,-0.4,0.5,-0.5), label = paste(species, "-", phylum))) +
  theme_light()

ggsave("figures/experimental_data.png", width = 12, height = 8)

write.csv(gb_selected_f, "data/species_lab_limits.csv", row.names = F)
