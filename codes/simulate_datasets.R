source("codes/simulate_species.R")

# Set a seed
set.seed(2023)

# Create a folder to hold the datasets
fs::dir_create("data/sim_datasets")

# Set the number of datasets
n_datasets <- 1

# Set the number of species
n_species <- 5

# Set environmental parameters
mean_site <- 15
sd_site <- 10

# Set number of cells
n_cells <- 500

# Set the means for each species
sel_x_species <- c(29, 27, 26, 20, 15)

# And also the variance
sel_xhat_species <- c(2, 3, 1, 2, 4)

# Print the 95% niche for each combination
for (x in 1:n_species) cat((sel_x_species[x] - (2 * sel_xhat_species[x])), "-", (sel_x_species[x] + (2 * sel_xhat_species[x])), "\n")

# Simulate the datasets
for (ds in 1:n_datasets) {
  # Simulate the species
  dataset <- lapply(1:n_species, function(sp_i){
    simulated_sp <- sim_species(x_species = sel_x_species[sp_i],
                                xhat_species = sel_xhat_species[sp_i],
                                x_site = mean_site, 
                                xhat_site = sd_site,
                                ncells = n_cells)
    
    data.frame(species = paste0("species_", sp_i),
               sst = simulated_sp$surface,
               suitability = simulated_sp$surface_suitability,
               true_occurrence = simulated_sp$true_occurrence,
               occurrence = simulated_sp$sampled_occurrence)
  })
  dataset <- do.call("rbind", dataset)
  
  write.csv(dataset, paste0("data/sim_datasets/dataset_", ds, ".csv"),
            row.names = F)
}
