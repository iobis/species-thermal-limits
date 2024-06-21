# Load packages
library(rethinking)
library(dplyr)
library(ggplot2)

# Load one dataset
data_raw <- read.csv("data/sim_datasets/dataset_1.csv")

# Create a list to hold results
sim_out <- list()

# Set percentages of absences
priors <- expand.grid(
  prior_m = seq(10, 30, by = 5),
  prior_sd = seq(2, 5, by = 1)
)

# Make data object
dat <- list(
  N = nrow(data_raw),
  N_spp = length(unique(data_raw$species)),
  sid = as.integer(as.factor(data_raw$species)),
  sst = data_raw$sst,
  y = data_raw$occurrence,
  P_mu = 20,
  P_sd = 5
)

# Run
for (z in 1:nrow(priors)) {
  
  # Change priors
  dat$P_mu <- priors$prior_m[z]
  dat$P_sd <- priors$prior_sd[z]
  
  # max and sd for each species, from sim data
  spp_mu <- rep(NA,dat$N_spp)
  spp_sd <- spp_mu
  for ( i in 1:dat$N_spp ) {
    suit <- data_raw$suitability[ dat$sid==i ]
    sst <- data_raw$sst[ dat$sid==i ]
    spp_mu[i] <- sst[ suit==max(suit) ]
    sst_b <- data_raw$sst[dat$sid==i & dat$y == 1]
    spp_sd[i] <- sd(sst_b)
  }
  
  # run model
  
  m0 <- cstan( file="codes/model2.stan" , data=dat , rstan_out=FALSE )
  
  prec_res <- precis(m0,2)
  prec_res$expected <- c(spp_mu, spp_sd, 0)
  
  sim_out[[z]] <- prec_res
  
}

# Select only those with a certain SD prior
plots_data_f <- lapply(unique(priors$prior_sd), function(x){
  
  tg_sd <- x
  sim_out_sel <- sim_out[which(priors$prior_sd == tg_sd)]
  
  # Mean
  plots_data <- list()
  for (i in 1:5) {
    plots_data[[i]] <- lapply(sim_out_sel, function(x){
      x[i,c(1,3,4,7)]
    })
    plots_data[[i]] <- do.call("rbind", plots_data[[i]])
    plots_data[[i]]$mean_prior <- priors$prior_m[priors$prior_sd == tg_sd]
    plots_data[[i]]$species <- paste("Species", i)
    plots_data[[i]]$expected <- plots_data[[i]]$expected[1]
  }
  
  plots_data <- do.call("rbind", plots_data)
  plots_data$sd_prior <- tg_sd
  plots_data
})

plots_data_f <- do.call("rbind", plots_data_f)

ggplot(plots_data_f[plots_data_f$sd_prior == 5,]) +
  geom_line(aes(x = mean_prior, y = mean)) +
  geom_ribbon(aes(x = mean_prior, ymin = `5.5%`, ymax = `94.5%`), alpha = .2, fill = "#168aad") +
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


