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

# Identify species from GlobTHERM that are marine ----
gb <- read.csv("data-raw/GlobalTherm_upload_02_11_17.csv")

gb_species <- paste(gb$Genus, gb$Species)
gb_species <- trimws(iconv(gb_species, from = "ISO-8859-1", to = "UTF-8"),
                     whitespace = "[\\h\\v]")

gb_taxonomy <- obistools::match_taxa(gb_species)

sp_list <- cbind(species = gb_species, gb_taxonomy)
sp_list <- sp_list[!is.na(sp_list$scientificName),]


# Load OBIS full export to see number of records for each species
obis <- open_dataset("~/Research/mpa_europe/mpaeu_shared/obis_20231025.parquet")

obis_records <- obis %>%
  select(scientificName, phylum, class, order, family, AphiaID) %>%
  filter(AphiaID %in% sp_list$acceptedNameUsageID) %>%
  #filter(scientificName %in% sp_list$scientificName) %>%
  group_by(scientificName, phylum, class, order, family) %>%
  count() %>%
  collect()

# Check fish species
obis_records %>%
  arrange(desc(n)) %>%
  filter(class == "Teleostei") %>%
  View()
