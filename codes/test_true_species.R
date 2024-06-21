# Test with real species and a simulated version of the same species

# Load packages ----
library(rethinking)
library(terra)
source("codes/simulate_species.R")
set.seed(2023)

# Load data ----
eclu <- read.csv("data/echinometra_lucunter_pa.csv")
head(eclu)

# Simulate a species with similar conditions
# Load environmental layer
sst <- rast("data/env_data/BO21_tempmean_ss.tif")
summary(sst)
ncells <- nrow(as.data.frame(sst))

sim_eclu <- sim_species(
  x_species = mean(eclu$sst[eclu$presence == 1]),
  xhat_species = sd(eclu$sst[eclu$presence == 1]),
  x_site = global(sst, "mean", na.rm = T)[1,1],
  xhat_site = global(sst, "sd", na.rm = T)[1,1], 
  site_min = minmax(sst)[1,1],
  site_max = minmax(sst)[2,1],
  ncells = ncells,
  samp_prob = 0.5
)

# Get a similar number of presence/absence as the true one
presence <- which(sim_eclu$sampled_occurrence == 1)
absence <- which(sim_eclu$sampled_occurrence == 0)

n_presence <- sum(eclu$presence)
n_absence <- nrow(eclu) - n_presence

sim_eclu <- as.data.frame(sim_eclu)

sim_eclu_balanced <- sim_eclu[c(sample(presence, n_presence),
                                sample(absence, n_absence)),]
colnames(sim_eclu_balanced) <- c("sst", "suitability", "true_presence", "presence")

# Prepare true/simulated dataset
multi_data <- dplyr::bind_rows(
  cbind(eclu, species = "true_species"),
  cbind(sim_eclu_balanced, species = "simulated_species")
)


# Prepare data for stan model
dat <- list(
  N = nrow(multi_data),
  N_spp = length(unique(multi_data$species)),
  sid = as.integer(as.factor(multi_data$species)),
  sst = multi_data$sst,
  y = multi_data$presence
)


# run model
m0 <- cstan( file="codes/model1.stan" , data=dat , rstan_out=FALSE )

precis(m0,2) # Simulated is first

# max and sd for each species, from sim data
spp_mu <- rep(NA, dat$N_spp)

spp_mu[2] <- sim_eclu_balanced$sst[sim_eclu_balanced$suitability == max(sim_eclu_balanced$suitability)]

# True species mean from data
mean(eclu$sst[eclu$presence == 1])

# True species experimental optimum: 29.4
# True species experimental limit: 36

plot(density(eclu$sst[eclu$presence == 1]))

plot(y = dnorm(20:40, mean = 31.8, sd = 5), x = 20:40, type = "l")
abline(v = 36)
