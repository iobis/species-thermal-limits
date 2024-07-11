# Get temperature for fish species

# Load packages 
library(duckdb)
library(DBI)
library(terra)
library(h3r)
library(dplyr)
library(ggplot2)
library(patchwork)

# List selected species
edna_fishes <- read.csv("data/funcdiversity_edna_fishes.csv", sep = "\t")
names(edna_fishes)[1] <- "species"
edna_fishes_aphia <- read.csv("../../obis/research/marineheritage_sst/data/unique_species.csv")
edna_fishes <- left_join(edna_fishes, edna_fishes_aphia[,c("species", "AphiaID")])

sel_species <- na.omit(unique(edna_fishes$AphiaID))

# Set up duckdb connection and register sst table
sst <- sfarrow::st_read_parquet("data/temperature_grid.parquet")
sst <- sst %>%
  st_drop_geometry()

con <- dbConnect(duckdb())
dbSendQuery(con, "install httpfs; load httpfs;")
duckdb_register(con, "sst", sst)

# Join cells list and gridded species dataset
species <- dbGetQuery(con, paste0("
  select species, AphiaID, surface, bottom
  from sst
  inner join (
    select *
    from read_parquet('s3://obis-products/speciesgrids/h3_7/*')
    where AphiaID in (", paste(sel_species, collapse = ", "), ")
  ) h3 on sst.h3_7 = h3.h3_07
  group by species, AphiaID, surface, bottom
"))

DBI::dbDisconnect(con)


# Get counts and summarise
sp_count <- species %>%
  group_by(species) %>%
  count()

sp_count <- left_join(sp_count, edna_fishes[,c("species", "Habitat")])

sp_count <- sp_count %>%
  mutate(hab = case_when(Habitat %in% c("reef-associated", "pelagic", "pelagic-oceanic") ~ "surface",
                         .default = "bottom")) %>%
  # Limit by those with at least 10 cells
  filter(n >= 10)


surface_species <- species %>%
  filter(species %in% sp_count$species[sp_count$hab == "surface"]) %>%
  select(species, surface) %>%
  group_by(species) %>%
  summarise(mean = mean(surface),
            sd = sd(surface))

bottom_species <- species %>%
  filter(species %in% sp_count$species[sp_count$hab == "bottom"]) %>%
  select(species, bottom) %>%
  group_by(species) %>%
  summarise(mean = mean(bottom),
            sd = sd(bottom))

surface_species <- left_join(surface_species, sp_count[,c("species", "n")])
bottom_species <- left_join(bottom_species, sp_count[,c("species", "n")])

p1 <- ggplot(surface_species) +
  geom_point(aes(x = mean, y = sd, size = log(n)), alpha = .2, color = "#11b5ae") +
  geom_smooth(aes(x = mean, y = sd), color = "grey20") +
  ggtitle("Pelagic species") +
  theme_light()

p2 <- ggplot(bottom_species) +
  geom_point(aes(x = mean, y = sd, size = log(n)), alpha = .2, color = "#4046ca") +
  geom_smooth(aes(x = mean, y = sd), color = "grey20") +
  ggtitle("Benthic species") +
  theme_light()


p1 + p2
ggsave("figures/fishes_temperature.png", width = 12, height = 6)


















# Try to run the model
library(rethinking)

con <- dbConnect(duckdb())
dbSendQuery(con, "install httpfs; load httpfs;")
duckdb_register(con, "sst", sst)

sel_species <- na.omit(unique(edna_fishes$AphiaID[edna_fishes$species %in% surface_species$species]))

# Join cells list and gridded species dataset
species_h3 <- dbGetQuery(con, paste0("
  select species, AphiaID, surface, bottom, h3_07
  from sst
  inner join (
    select *
    from read_parquet('s3://obis-products/speciesgrids/h3_7/*')
    where AphiaID in (", paste(sel_species, collapse = ", "), ")
  ) h3 on sst.h3_7 = h3.h3_07
  group by species, AphiaID, surface, bottom, h3_07
"))

DBI::dbDisconnect(con)

# Prepare data for model
species_pa <- list()

for (i in 1:length(unique(species_h3$species))) {
  sp <- unique(species_h3$species)[i]
  
  sel_data <- species_h3[species_h3$species == sp,]
  other_data <- species_h3[species_h3$species != sp,]
  
  other_data <- other_data[!other_data$h3_07 %in% sel_data$h3_07,]
  
  near_h3 <- h3r::gridDisk(sel_data$h3_07, rep(10, length(sel_data$h3_07)))
  near_h3 <- unlist(near_h3)
  
  other_data <- other_data[other_data$h3_07 %in% near_h3,]
  
  species_pa[[i]] <- data.frame(rbind(
    cbind(sel_data[,c("species", "surface")], presence = 1),
    data.frame(species = sel_data$species[1],
               surface = other_data$surface,
               presence = 0)
  ))
  
}

species_pa <- bind_rows(species_pa)

dat <- list(
  N = nrow(species_pa),
  N_spp = length(unique(species_pa$species)),
  sid = as.integer(as.factor(species_pa$species)),
  sst = species_pa$surface,
  y = species_pa$presence
)

# run model
m0 <- cstan( file="codes/model1.stan" , data=dat , rstan_out=FALSE )

model_result <- precis(m0, 2)

msp <- levels(as.factor(species_pa$species))
surface_species[surface_species$species %in% msp,]











# Get SST data for unique H3
unique_h3 <- unique(species_data$h3_07)

# Get coordinates
species_coords <- h3r::cellToLatLng(unique_h3)

# Load SST data
sst <- rast("../mpaeu_sdm/data/env/current/thetao_baseline_depthsurf_mean.tif")
species_sst <- terra::extract(sst,
                              species_coords[,c(2,1)], # Important! Inverted lon/lat
                              ID = F, 
                              method = "bilinear")

# create a function to imput from nearby
get_near <- function(lon_lat, layer) {
  
  cells <- cellFromXY(layer, lon_lat)
  
  near_cells <- terra::adjacent(layer, cells, directions = "16")
  
  near_data <- apply(near_cells, 1, function(x){
    values <- terra::extract(layer, x)
    mean(values[,1], na.rm = T)
  })
  
  return(near_data)
}

# Input nearby
sum(is.na(species_sst[,1]))
species_sst[is.na(species_sst[,1]), 1] <- get_near(species_coords[is.na(species_sst[,1]), c(2,1)],
                                                   sst)

sum(is.na(species_sst[,1])) # Fair reduction...

# Join with species data
colnames(species_sst) <- "SST"
species_sst <- cbind(species_sst, h3_07 = unique_h3)
species_data <- left_join(species_data, species_sst)

# Summarise
quantile_df <- function(x) {
  tibble(
    metric = c("mean", "sd", "q_0.05", "q_0.5", "q_0.95"),
    value = c(
      mean(x, na.rm = T), sd(x, na.rm = T),
      quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE)
    )
  )
}

thermal_limits <- species_data %>%
  group_by(species, AphiaID) %>%
  reframe(quantile_df(SST))

head(thermal_limits)