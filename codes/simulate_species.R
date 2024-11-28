#' Simulate a virtual species distribution
#'
#' @param x_species mean of the parameter for the species
#' @param xhat_species SD of the parameter for the species
#' @param x_site mean of the parameter for the site (sampling surface)
#' @param xhat_site SD of the parameter for the site (sampling surface)
#' @param site_min the minimum value that the site can have (use `-Inf` if you want to avoid truncation)
#' @param site_max the maximum value that the site can have (use `Inf` if you want to avoid truncation)
#' @param ncells number of cells available for the species
#' @param samp_prob probability of sampling/finding a species where it occurs
#' @param samp_equal if `TRUE`, the final values (i.e. samples) will be sampled to ensure equal number
#'   of presences (1) and absences (0). Ensure a large number of ncells, so that there is enough points to
#'   sample. Otherwise, the function will return an error.
#' @param samp_equal_target the number of presences and absences that should be sampled when `samp_equal = T`
#'
#' @return a list with sampled and true occurrence
#' @export
#'
#' @examples
#' sim_results <- sim_species()
#'
#' Our hypothetical species is tropical, with its optimum on 29 degrees.
#' We add a SD of 2, such that 95% of the species niche are between
#' 25 degrees and 33 degrees (a relatively narrow niche)

#' sp_mean <- 29
#' sp_sd <- 2

#' For the probability of sampling the species we are leaving the default = 0.5
#' We assume that the species is never found in places where it does not occur
#' (i.e, we are not assuming any false positives)

#' sim_results <- sim_species(x_species = sp_mean,
#'                            xhat_species = sp_sd,
#'                            x_site = 15, xhat_site = 10)

#' Table true occurrence vs sampled

#' table(true = sim_results$true_occurrence,
#'       sampled = sim_results$sampled_occurrence)

#' Plot

#' plot(x = sim_results$surface, y = sim_results$surface_suitability,
#'      col = c("#3B7B9A55", "#EFCC3055")[as.factor(sim_results$true_occurrence)],
#'      pch = 20, cex = 2, xlab = "SST", ylab = "Suitability",
#'      xlim = c(3, 35))
#' points(x = sim_results$surface, y = sim_results$surface_suitability,
#'      col = c("#3B7B9A55", "#EFCC3055")[as.factor(sim_results$sampled_occurrence)],
#'      pch = 1, cex = 3)
#' abline(v = max(sim_results$surface), col = "black")
#' abline(v = sp_mean, lty = 2, col = "black")
#' abline(v = c(sp_mean - 2*sp_sd, sp_mean + 2*sp_sd), lty = 2, col = "grey80")
#' 
sim_species <- function(x_species = 29, xhat_species = 2,
                        x_site = 28, xhat_site = 8,
                        site_min = 0, site_max = 30,
                        ncells = 100, samp_prob = 0.5,
                        samp_equal = TRUE, samp_equal_target = 20){
  
  # We start by creating the surface, i.e. the sites. We assume that the distribution
  # of the sites follow a Gaussian distribution.
  # We want to have a truncated distribution. Thus the available temperatures
  # will go up to a certain value, while the species can in fact survive temperatures
  # higher than this threshold.
  samp_surf_const <- truncnorm::rtruncnorm(
    n = ncells, a = site_min, b = site_max, mean = x_site, sd = xhat_site
  )
  
  # We convert the surface to suitability values based again on a Gaussian distribution
  # with mean = x_species and sd = xhat_species
  suitability <- dnorm(samp_surf_const, mean = x_species, sd = xhat_species)
  
  # Normalize to 0-1 for better plotting
  suitability <- (suitability - min(suitability)) / (max(suitability) - min(suitability))
  
  # Convert to presence/absence of the species based on binomial distribution, with
  # probability equal to the suitability.
  true_occurrence <- rbinom(length(suitability), size = 1, prob = suitability)
  
  # Sample the occurrence whith binomial probability, simulating a sampling campaign
  sampled_occurrence <- rbinom(length(true_occurrence), size = 1, prob = true_occurrence * samp_prob)

  # If equal numbers
  if (samp_equal) {
    sampled_occurrence <- as.integer(sampled_occurrence)
    freqs <- as.data.frame(table(sampled_occurrence))
    if (nrow(freqs) == 1) stop("Only one class was generated.")
    if (any(freqs$Freq < samp_equal_target)) {
      stop("Not enough number of points in one or more classes. Increasing the number of cells might help.")
    }
    sampled_occurrence_sample <- tapply(seq_len(length(sampled_occurrence)),
                                        as.factor(sampled_occurrence),
                                        sample, size = samp_equal_target)
    sampled_occurrence_sample <- unname(unlist(sampled_occurrence_sample))

    true_occurrence <- true_occurrence[sampled_occurrence_sample]
    sampled_occurrence <- sampled_occurrence[sampled_occurrence_sample]
    suitability <- suitability[sampled_occurrence_sample]
    samp_surf_const <- samp_surf_const[sampled_occurrence_sample]

  }
  
  return(list(
    surface = samp_surf_const,
    surface_suitability = suitability,
    true_occurrence = as.integer(true_occurrence),
    sampled_occurrence = as.integer(sampled_occurrence)
  ))
  
}


#' Simulate a virtual species distribution based on true environmental data
#'
#' @param x_species mean of the parameter for the species
#' @param xhat_species SD of the parameter for the species
#' @param site_values true values of environmental data 
#' @param site_min optional, the minimum value that the site can have, to truncate the true values
#' @param site_max optional, the maximum value that the site can have, to truncate the true values
#' @param ncells optional, number of cells available for the species. If more values are supplied
#'  in `site_values`, then those are sampled to `ncells`
#' @param samp_prob probability of sampling/finding a species where it occurs
#'
#' @return a list with sampled and true occurrence
#' @export
#'
#' @examples
#' sim_results <- sim_species_env()
#' 
sim_species_env <- function(x_species = 29, xhat_species = 2,
                            site_values,
                            site_min = NULL, site_max = NULL,
                            ncells = NULL, samp_prob = 0.5){
  
  site_values <- na.omit(site_values)

  if (!is.null(site_min)) {
    site_values <- site_values[site_values >= site_min]
  }
  if (!is.null(site_max)) {
    site_values <- site_values[site_values <= site_max]
  }
  if (!is.null(ncells) && length(site_values) > ncells) {
    site_values <- site_values[site_values <= site_max]
  } else if (!is.null(ncells)) {
    cat("`site_values` is",
     ifelse(length(site_values) < ncells, "less than `ncells`.", "equal than `ncells`."),
     "Ignoring `ncells`.\n")
  }
  
  # We convert the surface to suitability values based again on a Gaussian distribution
  # with mean = x_species and sd = xhat_species
  suitability <- dnorm(site_values, mean = x_species, sd = xhat_species)
  
  # Normalize to 0-1 for better plotting
  suitability <- (suitability - min(suitability)) / (max(suitability) - min(suitability))
  
  # Convert to presence/absence of the species based on binomial distribution, with
  # probability equal to the suitability.
  true_occurrence <- rbinom(length(suitability), size = 1, prob = suitability)
  
  # Sample the occurrence whith binomial probability, simulating a sampling campaign
  sampled_occurrence <- rbinom(length(true_occurrence), size = 1, prob = true_occurrence * samp_prob)
  
  return(list(
    surface = site_values,
    surface_suitability = suitability,
    true_occurrence = as.integer(true_occurrence),
    sampled_occurrence = as.integer(sampled_occurrence)
  ))
  
}