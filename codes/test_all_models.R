library(rethinking)
library(ggplot2)
source("functions/sim_dataset.R")
set.seed(20)

n_species <- 10
n_pts <- 60

acanthurus <- rfishbase::species()
acanthurus <- acanthurus[acanthurus$Genus == "Acanthurus",]

acnth_estimates <- rfishbase::estimate()
acnth_estimates <- acnth_estimates[acnth_estimates$SpecCode %in% acanthurus$SpecCode,]

acnth_pref <- na.omit(acnth_estimates$TempPrefMean)
acnth_pref_m <- median(acnth_pref)
acnth_pref_sd <- sd(acnth_pref)

opt_t <- acnth_pref_m
max_cov <- 0.1

dataset_full <- sim_dataset_phylo(
    N = n_pts * n_species,                       # Number of surveys
    n_sp = n_species,                            # Number of species
    p = 0.8, 
    global_optimum = log(opt_t),
    max_cov = max_cov,
    rho_val = 1,
    equal_n_pres = TRUE
)
aggregate(dataset_full$data$y, list(dataset_full$data$sid), table)

# Generate a dataset with some species with low records
dataset_reduced <- dataset_full
to_reduce <- sample(unique(dataset_reduced$data$sid), 4)
n_recs <- sample(1:5, length(to_reduce), replace = T)
to_remove <- vector("list", length(to_reduce))
for (i in seq_along(to_reduce)) {
    ab <- which(dataset_reduced$data$sid == to_reduce[i] & dataset_reduced$data$y == 0)
    pr <- which(dataset_reduced$data$sid == to_reduce[i] & dataset_reduced$data$y == 1)
    to_remove[[i]] <- c(ab, sample(pr, (n_pts/2) - n_recs[i]))
}
to_remove <- do.call("c", to_remove)
dataset_reduced$data$sid <- dataset_reduced$data$sid[-to_remove]
dataset_reduced$data$y <- dataset_reduced$data$y[-to_remove]
dataset_reduced$data$sst <- dataset_reduced$data$sst[-to_remove]
dataset_reduced$data$N <- length(dataset_reduced$data$sid)

# Test model 1
m1_full <- cstan(file = paste0("codes/model1.stan"),
                 data = dataset_full$data[c("N", "N_spp", "sid", "sst", "y")],
                 rstan_out = FALSE)

m1_reduced <- cstan(file = paste0("codes/model1.stan"),
                    data = dataset_reduced$data[c("N", "N_spp", "sid", "sst", "y")],
                    rstan_out = FALSE)

# Test model 4
m4_full <- cstan(file = paste0("codes/model4.stan"),
                 data = dataset_full$data[c("N", "N_spp", "sid", "sst", "y")],
                 rstan_out = FALSE)

m4_reduced <- cstan(file = paste0("codes/model4.stan"),
                    data = dataset_reduced$data[c("N", "N_spp", "sid", "sst", "y")],
                    rstan_out = FALSE)

# Test model 5
m5_full <- cstan(file = paste0("codes/model5.stan"),
                 data = dataset_full$data[c("N", "N_spp", "sid", "sst", "y")],
                 rstan_out = FALSE)

m5_reduced <- cstan(file = paste0("codes/model5.stan"),
                    data = dataset_reduced$data[c("N", "N_spp", "sid", "sst", "y")],
                    rstan_out = FALSE)

# Test model 6
m6_full <- cstan(file = paste0("codes/model6.stan"),
                 data = dataset_full$data[c("N", "N_spp", "sid", "sst", "y")],
                 rstan_out = FALSE)

m6_reduced <- cstan(file = paste0("codes/model6.stan"),
                    data = dataset_reduced$data[c("N", "N_spp", "sid", "sst", "y")],
                    rstan_out = FALSE)

# Test model 7
m7_full <- cstan(file = paste0("codes/model7.stan"),
                 data = dataset_full$data[c("N", "N_spp", "sid", "sst", "y")],
                 rstan_out = FALSE)

m7_reduced <- cstan(file = paste0("codes/model7.stan"),
                    data = dataset_reduced$data[c("N", "N_spp", "sid", "sst", "y")],
                    rstan_out = FALSE)

# Test model 8
m8_full <- cstan(file = paste0("codes/model8.stan"),
                 data = dataset_full$data,
                 rstan_out = FALSE)

m8_reduced <- cstan(file = paste0("codes/model8.stan"),
                    data = dataset_reduced$data,
                    rstan_out = FALSE)

# Plot
results <- lapply(c(1, 4, 5, 7, 8), \(id) {
    df <- data.frame(model = id, sid = 1:n_species, mean = NA, q25 = NA, q75 = NA)
    if (id < 8) {
        mr <- eval(parse(text = paste0("mr <- precis(m", id, "_full, 2, prob = .50)")))
        mr <- mr[grepl("tmu\\[", rownames(mr)), ]
        if (any(grepl("log", rownames(mr)))) {
            mr <- exp(mr)
        }
        df$mean <- mr$mean
        df$q25 <- mr$`25%`
        df$q75 <- mr$`75%`
    } else {
        mr <- eval(parse(text = paste0("mr <- extract.samples(m", id, "_full)")))
        true_values <- exp(dataset_full$raw$f[, 1])
        estimated_values <- data.frame(
            q_025 = exp(apply(mr$f, 2:3, quantile, .25))[,1],
            mean = exp(apply(mr$f, 2:3, mean)[,1]),
            q_075 = exp(apply(mr$f, 2:3, quantile, .75))[,1]
        )
        df$mean <- estimated_values$mean
        df$q25 <- estimated_values$q_025
        df$q75 <- estimated_values$q_075
    }
    return(df)
})

results <- do.call("rbind", results)
true_sst <- exp(dataset_full$raw$f[, 1])
results$expected <- rep(true_sst, length(unique(results$model)))

ggplot(results) +
    geom_pointrange(
        aes(x = expected, y = mean, ymin = q25, ymax = q75),
        color = "#1860a8"
    ) +
    geom_abline(color = "grey80", linetype = 2) +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~model) +
    xlab("Expected") + ylab("Estimated") + 
    ggtitle("Full data",
            glue::glue("N pres = {n_pts/2}, N abs = {n_pts/2}, N species = {n_species}"))


results_red <- lapply(c(1, 4, 5, 7, 8), \(id) {
    df <- data.frame(model = id, sid = 1:n_species, mean = NA, q25 = NA, q75 = NA)
    if (id < 8) {
        mr <- eval(parse(text = paste0("mr <- precis(m", id, "_reduced, 2, prob = .50)")))
        mr <- mr[grepl("tmu\\[", rownames(mr)), ]
        if (any(grepl("log", rownames(mr)))) {
            mr <- exp(mr)
        }
        df$mean <- mr$mean
        df$q25 <- mr$`25%`
        df$q75 <- mr$`75%`
    } else {
        mr <- eval(parse(text = paste0("mr <- extract.samples(m", id, "_reduced)")))
        true_values <- exp(dataset_reduced$raw$f[, 1])
        estimated_values <- data.frame(
            q_025 = exp(apply(mr$f, 2:3, quantile, .25))[,1],
            mean = exp(apply(mr$f, 2:3, mean)[,1]),
            q_075 = exp(apply(mr$f, 2:3, quantile, .75))[,1]
        )
        df$mean <- estimated_values$mean
        df$q25 <- estimated_values$q_025
        df$q75 <- estimated_values$q_075
    }
    return(df)
})

results_red <- do.call("rbind", results_red)
true_sst <- exp(dataset_reduced$raw$f[, 1])
results_red$expected <- rep(true_sst, length(unique(results_red$model)))
n_pts_ind <- aggregate(dataset_reduced$data$y, list(dataset_reduced$data$sid), length)
results_red$n_pts <- rep(n_pts_ind$x, length(unique(results_red$model)))

ggplot(results_red) +
    geom_pointrange(
        aes(x = expected, y = mean, ymin = q25, ymax = q75),
        color = "#1860a8"
    ) +
    geom_text(
        aes(x = expected, y = q75 + .8, label = n_pts), size = rel(3)
    ) +
    geom_abline(color = "grey80", linetype = 2) +
    theme_light() +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~model) +
    xlab("Expected") + ylab("Estimated") + 
    ggtitle("Reduced data",
            glue::glue("N pres = {n_pts/2}, N abs = {n_pts/2}, N species = {n_species}"))
