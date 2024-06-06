# Load packages
library(rethinking)
library(dplyr)
library(ggplot2)

# Load one dataset
data_raw <- read.csv("data/sim_datasets/dataset_1.csv")

# Create a list to hold results
sim_out <- list()

# Set percentages of absences
percentages <- seq(0, 1, by = 0.1)

# Run
for (z in 1:length(percentages)) {
  
  # Separate occurrences
  data_occurrences <- data_raw %>%
    filter(occurrence == 1)
  
  # Separate absences
  data_absences <- data_raw %>%
    filter(occurrence == 0) %>%
    group_by(species) %>%
    slice_sample(prop = percentages[z])
  
  #print(data_absences %>% group_by(species) %>% count())
  # Bind rows
  data_edited <- rbind(data_occurrences, data_absences)
  
  # Make data object
  dat <- list(
    N = nrow(data_edited),
    N_spp = length(unique(data_edited$species)),
    sid = as.integer(as.factor(data_edited$species)),
    sst = data_edited$sst,
    y = data_edited$occurrence
  )
  
  # max and sd for each species, from sim data
  spp_mu <- rep(NA,dat$N_spp)
  spp_sd <- spp_mu
  for ( i in 1:dat$N_spp ) {
    suit <- data_edited$suitability[ dat$sid==i ]
    sst <- data_edited$sst[ dat$sid==i ]
    spp_mu[i] <- sst[ suit==max(suit) ]
    sst_b <- data_edited$sst[dat$sid==i & dat$y == 1]
    spp_sd[i] <- sd(sst_b)
  }
  
  # run model
  
  m0 <- cstan( file="codes/model1.stan" , data=dat , rstan_out=FALSE )
  
  prec_res <- precis(m0,2)
  prec_res$expected <- c(spp_mu, spp_sd, 0)
  
  sim_out[[z]] <- prec_res
  
}

# Mean
plots_data <- list()
for (i in 1:5) {
  plots_data[[i]] <- lapply(sim_out, function(x){
    x[i,c(1,3,4,7)]
  })
  plots_data[[i]] <- do.call("rbind", plots_data[[i]])
  plots_data[[i]]$percentages <- percentages
  plots_data[[i]]$species <- paste("Species", i)
  plots_data[[i]]$expected <- plots_data[[i]]$expected[1]
}

plots_data <- do.call("rbind", plots_data)

ggplot(plots_data) +
  geom_line(aes(x = percentages, y = mean)) +
  geom_ribbon(aes(x = percentages, ymin = `5.5%`, ymax = `94.5%`), alpha = .2, fill = "#168aad") +
  geom_hline(aes(yintercept = expected), linetype = 2) +
  facet_wrap(~ species, scales = "free_y") +
  theme_light()



# SD
plots_data_sd <- list()
for (i in 6:10) {
  plots_data_sd[[i]] <- lapply(sim_out, function(x){
    x[i,c(1,3,4,7)]
  })
  plots_data_sd[[i]] <- do.call("rbind", plots_data_sd[[i]])
  plots_data_sd[[i]]$percentages <- percentages
  plots_data_sd[[i]]$species <- paste("Species", (i - 5))
  plots_data_sd[[i]]$expected <- plots_data_sd[[i]]$expected
}

plots_data_sd <- do.call("rbind", plots_data_sd)

ggplot(plots_data_sd) +
  geom_line(aes(x = percentages, y = mean)) +
  geom_ribbon(aes(x = percentages, ymin = `5.5%`, ymax = `94.5%`), alpha = .2, fill = "#168aad") +
  geom_hline(aes(yintercept = expected), linetype = 2) +
  facet_wrap(~ species, scales = "free_y") +
  theme_light()


