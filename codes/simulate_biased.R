set.seed(2024)
source("codes/simulate_species.R")

# Our hypothetical species is tropical, with its optimum on 28
sp_mean <- 28
sp_sd <- 3

# We start with the normal simulation
sim_results <- sim_species(x_species = sp_mean,
                           xhat_species = sp_sd,
                           x_site = 26, xhat_site = 10, ncells = 500)

# We have a surface of temperatures, with its suitability and true presence/absence
# We will now simulate sampling campaigns, biased or not towards a certain temperature
# As done in the simulation, we will use an error probability of 0.5
# That is, sampling a presence as an absence.
prob_error <- 0.5

# For the biased, we will push the sampling towards a certain temperature
bias_factor <- 0.25
collection_prob <- dnorm(sim_results$surface, mean = 22, sd = 2)
collection_prob <- bias_factor + (1-bias_factor) * (collection_prob - min(collection_prob)) / (max(collection_prob) - min(collection_prob))
# Sites with temperatures around 20 will have up to 75% more chance of being sampled

# Now simulate the sampling
total <- 200
normal_sampling <- sample(1:length(sim_results$true_occurrence), total, replace = F)

biased_sampling <- sample(1:length(sim_results$true_occurrence), total, prob = collection_prob, replace = F)


# Get the values
data <- data.frame(do.call("cbind", sim_results))[,1:3]

normal_data <- data[normal_sampling,]
biased_data <- data[biased_sampling,]

# Add the "sampled" occurrence (with possible error)
normal_data$sampled_occurrence <- rbinom(length(normal_data$true_occurrence),
                                         size = 1, prob = normal_data$true_occurrence * prob_error)

table(true = normal_data$true_occurrence,
      sampled = normal_data$sampled_occurrence)


biased_data$sampled_occurrence <- rbinom(length(biased_data$true_occurrence),
                                         size = 1, prob = biased_data$true_occurrence * prob_error)

table(true = biased_data$true_occurrence,
      sampled = biased_data$sampled_occurrence)

# Save datasets
normal_data$type <- "normal"
biased_data$type <- "biased"

write.csv(rbind(normal_data, biased_data), "data/sim_datasets/biased_dataset1.csv", row.names = F)


# Plot
plot(x = sim_results$surface, y = sim_results$surface_suitability,
     type = "l", col = "#3B7B9A", xlab = "SST", ylab = "Suitability", xlim = c(3, 35))
lines(x = sim_results$surface, y = collection_prob, type = "l", col = "#EFCC30")
abline(v = c(sp_mean - 2*sp_sd, sp_mean + 2*sp_sd), lty = 2, col = "grey80")
abline(v = sp_mean, lty = 2, col = "black")
text(x = sp_mean+0.5, y = 0.5, label = "true\nmean", col = "grey20", adj = 0)

samp_mean <- mean(normal_data$surface[normal_data$sampled_occurrence == 1])
abline(v = samp_mean, col = "black")
text(x = samp_mean+0.5, y = 0.3, label = "normal sampling\nmean", col = "grey20", adj = 0)

bias_mean <- mean(biased_data$surface[biased_data$sampled_occurrence == 1])
abline(v = mean(bias_mean), col = "#480ca8", lwd = 2)
text(x = mean(bias_mean)+0.5, y = 0.1, label = "biased sampling\nmean", col = "#480ca8", adj = 0)

legend("topleft", legend = c("Suitability", "Biased sampling probability"), col = c("#3B7B9A", "#EFCC30"), lty = 1)

