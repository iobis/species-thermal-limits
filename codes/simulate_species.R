#' Simulate a virtual species distribution
#'
#' @param x_species mean of the parameter for the species
#' @param xhat_species SD of the parameter for the species
#' @param x_site mean of the parameter for the site (sampling surface)
#' @param xhat_site SD of the parameter for the site (sampling surface)
#' @param site_min the minimum value that the site can have (use a very low value if you want to avoid truncation)
#' @param site_max the maximum value that the site can have (use a very high value if you want to avoid truncation)
#' @param ncells number of cells available for the species
#' @param samp_prob probability of sampling/finding a species where it occurs
#'
#' @return a list with sampled and true occurrence
#' @export
#'
#' @examples
#' sim_results <- sim_species()
#' 
sim_species <- function(x_species = 29, xhat_species = 2,
                        x_site = 28, xhat_site = 8,
                        site_min = 0, site_max = 30,
                        ncells = 100, samp_prob = 0.5){
  
  # We start by creating the surface, i.e. the sites. We assume that the distribution
  # of the sites follow a Gaussian distribution.
  # Because we will truncate this distribution, we produce more cells than requested
  samp_surf <- rnorm(ncells*2, mean = x_site, sd = xhat_site)
  
  # We want to have a truncated distribution. Thus the available temperatures
  # will go up to a certain value, while the species can in fact survive temperatures
  # higher than this threshold.
  samp_surf_const <- samp_surf[samp_surf > site_min & samp_surf < site_max]
  
  # We then sample to have the number of cells expected
  samp_surf_sel <- sample(samp_surf_const, ncells)
  samp_surf_sel <- samp_surf_sel[order(samp_surf_sel)]
  
  # We convert the surface to suitability values based again on a Gaussian distribution
  # with mean = x_species and sd = xhat_species
  suitability <- dnorm(samp_surf_sel, mean = x_species, sd = xhat_species)
  
  # Normalize to 0-1 for better plotting
  suitability <- (suitability - min(suitability)) / (max(suitability) - min(suitability))
  
  # Convert to presence/absence of the species based on binomial distribution, with
  # probability equal to the suitability.
  true_occurrence <- rbinom(length(suitability), size = 1, prob = suitability)
  
  # Sample the occurrence whith binomial probability, simulating a sampling campaign
  sampled_occurrence <- rbinom(length(true_occurrence), size = 1, prob = true_occurrence * samp_prob)
  
  return(list(
    surface = samp_surf_sel,
    surface_suitability = suitability,
    true_occurrence = true_occurrence,
    sampled_occurrence = sampled_occurrence
  ))
  
}


# Our hypothetical species is tropical, with its optimum on 29 degrees.
# We add a SD of 2, such that 95% of the species niche are between
# 25 degrees and 33 degrees (a relatively narrow niche)
sp_mean <- 29
sp_sd <- 2

# For the probability of sampling the species we are leaving the default = 0.5
# We assume that the species is never found in places where it does not occur
# (i.e, we are not assuming any false positives)
sim_results <- sim_species(x_species = sp_mean,
                           xhat_species = sp_sd,
                           x_site = 15, xhat_site = 10)

# Table true occurrence vs sampled
table(true = sim_results$true_occurrence,
      sampled = sim_results$sampled_occurrence)

# Plot
plot(x = sim_results$surface, y = sim_results$surface_suitability,
     col = c("#3B7B9A55", "#EFCC3055")[as.factor(sim_results$true_occurrence)],
     pch = 20, cex = 2, xlab = "SST", ylab = "Suitability",
     xlim = c(3, 35))
points(x = sim_results$surface, y = sim_results$surface_suitability,
     col = c("#3B7B9A55", "#EFCC3055")[as.factor(sim_results$sampled_occurrence)],
     pch = 1, cex = 3)
abline(v = max(sim_results$surface), col = "black")
abline(v = sp_mean, lty = 2, col = "black")
abline(v = c(sp_mean - 2*sp_sd, sp_mean + 2*sp_sd), lty = 2, col = "grey80")
