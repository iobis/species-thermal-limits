# prep data for Stan

library(rethinking)

data_raw <- read.csv("data/sim_datasets/dataset_1.csv")

dat <- list(
N = nrow(data_raw),
N_spp = length(unique(data_raw$species)),
sid = as.integer(as.factor(data_raw$species)),
sst = data_raw$sst,
y = data_raw$occurrence
)

# max and sd for each species, from sim data
spp_mu <- rep(NA,dat$N_spp)
spp_sd <- spp_mu  
for ( i in 1:dat$N_spp ) {
    suit <- data_raw$suitability[ dat$sid==i ]
    sst <- data_raw$sst[ dat$sid==i ]
    spp_mu[i] <- sst[ suit==max(suit) ]
}

# run model

m0 <- cstan( file="codes/model1.stan" , data=dat , rstan_out=FALSE )

precis(m0,2)

spp_mu

