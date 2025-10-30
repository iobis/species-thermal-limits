########### Thermal limits of marine species based on occurrence data ##########
# July of 2025
# Authors: Silas Principe, Richard McElreath, Pieter Provoost
# Contact: s.principe@unesco.org
#
############################### Tests for model 1 ##############################

# We start by creating a function which will enable us to simulate the data
sim_dataset <- function(
    N = 500,                                    # Number of surveys
    n_sp = 5,                                   # Number of species
    p = rbeta(1, 2, 2),                         # Detection probability
    mu_tmu = 20,                                # Global mean for tmu
    sigma_tmu = 3,                              # Global sd for tmu
    mu_tsd = 5,                                 # Global mean for tsd
    sigma_tsd = 1,                              # Global sd for tsd 
    tomax = rbeta(n_sp, 5, 1),                  # Max occupancy prob per species - you can also pass a fixed value
    site_min = -Inf,                            # Minimum of site (for truncation)
    site_max = Inf,                             # Maximum of site (for truncation)
    n_presence = NULL,                          # Number of presences
    n_absence = NULL,                           # Number of absences
    perc_sst_inc = 0                            # Percentage of increase in SST SD
) {

    tmu <- round(rnorm(n_sp, mu_tmu, sigma_tmu), 1)
    tsd <- round(abs(rnorm(n_sp, mu_tsd, sigma_tsd)), 1)

    if (length(tomax) == 1) {
        tomax <- rep(tomax, n_sp)
    }

    if (!is.null(n_presence) || !is.null(n_absence)) {
        if (is.null(n_presence)) {
            n_presence <- 10
        }
        if (is.null(n_absence)) {
            n_absence <- 10
        }
        max_try <- 1000
        species_dataset <- vector("list", n_sp)
        for (sp in seq_len(n_sp)) {
            condition_met <- FALSE
            N_try <- (n_presence + n_absence) * 10
            pres_list <- abs_list <- vector(mode = "list", length = max_try)
            k <- 1
            while(!condition_met & k <= max_try) {
                dataset <- .run_sim_ds(N = N_try, n_sp = 1, tmu[sp], tsd[sp],
                                       tomax[sp], p, site_min, site_max, perc_sst_inc)
                pres_list[[k]] <- dataset[dataset$y == 1,]
                abs_list[[k]] <- dataset[dataset$y == 0,]
                total_pres <- nrow(do.call("rbind", pres_list))
                total_abs <- nrow(do.call("rbind", abs_list))
                if (total_pres >= n_presence & total_abs >= n_absence) {
                    condition_met <- TRUE
                    pres_list <- do.call("rbind", pres_list)
                    abs_list <- do.call("rbind", abs_list)
                    if (total_pres > n_presence) {
                        pres_list <- pres_list[sample(seq_len(nrow(pres_list)), n_presence),]
                    }
                    if (total_abs > n_absence) {
                        abs_list <- abs_list[sample(seq_len(nrow(abs_list)), n_presence),]
                    }
                }
            }
            if (!condition_met) stop("Failed to produce dataset with this number of presence/absence")
            species_dataset[[sp]] <- rbind(abs_list, pres_list)
            species_dataset[[sp]]$sid <- sp
        }
        dataset <- as.data.frame(do.call("rbind", species_dataset))
        rownames(dataset) <- seq_len(nrow(dataset))
    } else {
        dataset <- .run_sim_ds(N, n_sp, tmu, tsd, tomax, p, site_min, site_max, perc_sst_inc)
    }

    list(
        dataset = dataset,
        p = p,
        tmu = tmu,
        tsd = tsd,
        tomax = tomax,
        mu_tmu = mu_tmu,
        sigma_tmu = sigma_tmu,
        mu_tsd = mu_tsd,
        sigma_tsd = sigma_tsd
    )
}

.run_sim_ds <- function(N, n_sp, tmu, tsd, tomax, p, site_min, site_max, perc_sst_inc) {
    # Species IDs
    sid <- rep(seq_len(n_sp), each = N / n_sp)

    # Simulate SST for each survey
    if (site_min != -Inf || site_max != Inf) {
        sst <- truncnorm::rtruncnorm(n = N, a = site_min, b = site_max, mean = tmu[sid], sd = tsd[sid])
    } else {
        sst <- rnorm(N, mean = tmu[sid], sd = (tsd[sid] + (tsd[sid] * perc_sst_inc)))
    }

    # Calculate occupancy probability (Gaussian suitability * tomax)
    q <- exp(-0.5 * ((sst - tmu[sid]) / tsd[sid])^2) * tomax[sid]

    # Simulate detections
    y <- numeric(N)
    occupancy <- numeric(N)
    for (i in seq_len(N)) {
        prob_occupied <- q[i] # Probability of occupancy given suitability
        occupancy[i] <- prob_occupied
        # Detected if occupied
        # prob_detected <- prob_occupied #presence is always true, otherwise use next line
        # prob_detected <- prob_occupied * p
        # Absence or undetected
        prob_absent <- (1 - prob_occupied) + prob_occupied * (1 - p)
        y[i] <- sample(c(1, 0), size = 1, prob = c(prob_occupied * p, prob_absent))
    }

    # Create dataset
    dataset <- data.frame(sid = sid, sst = sst, y = y, prob_occ = occupancy)

    return(dataset)
}


get_dataset_stats <- function(dataset) {
    tmu <- dataset$tmu
    tsd <- dataset$tsd
    dataset <- dataset$dataset
    agg <- tapply(dataset, dataset$sid, \(x) {
        data.frame(
            average_temperature = round(mean(x$sst), 1),
            sd_temperature = round(sd(x$sst), 1),
            number_presences = sum(x$y),
            number_absences = length(x$y) - sum(x$y)
        )
    })
    agg <- do.call("rbind", agg)
    agg$tmu <- tmu
    agg$tsd <- tsd
    rownames(agg) <- paste("Species", unique(dataset$sid))
    return(agg)
}

prepare_data_stan <- function(dataset) {
    dat <- list(
        N = nrow(dataset$dataset),
        n_sp = length(unique(dataset$dataset$sid)),
        sid = as.integer(as.factor(dataset$dataset$sid)),
        sst = dataset$dataset$sst,
        y = dataset$dataset$y
    )
    return(dat)
}

extract_precis <- function(model, dataset, transform = FALSE, change_names = FALSE) {
    precis_result <- precis(model, 2)
    rnames <- rownames(precis_result)
    #rownames(precis_result) <- NULL
    precis_result$what <- gsub("\\[.*.\\]", "", rnames)
    precis_result$sid <- NA
    precis_result$sid[grepl("tmu", precis_result$what)] <- unique(dataset$dataset$sid)
    precis_result$sid[grepl("tsd", precis_result$what)] <- unique(dataset$dataset$sid)
    precis_result$sid[grepl("tomax", precis_result$what)] <- unique(dataset$dataset$sid)
    precis_result$expected <- NA
    precis_result$expected[grepl("tmu", precis_result$what)] <- dataset$tmu
    precis_result$expected[grepl("tsd", precis_result$what)] <- dataset$tsd
    if (!is.null(dataset$p)) {
        precis_result$expected[grepl("p$", precis_result$what)] <- dataset$p
    }
    if (!is.null(dataset$tomax) & any(grepl("tomax", precis_result$what))) {
        precis_result$expected[grepl("tomax", precis_result$what)] <- dataset$tomax
    }
    if (any(grepl("spp_", rnames))) {
        if (any(grepl("sigma", rnames))) {
            precis_result$expected[grepl("spp_", precis_result$what)] <- c(
                dataset$mu_tmu, dataset$mu_tsd, NA
            )
            precis_result$expected[grepl("sigma", precis_result$what)] <- c(
                dataset$sigma_tmu, dataset$sigma_tsd, NA
            )
        } else {
            precis_result$expected[grepl("spp_", precis_result$what)] <- c(
                dataset$mu_tmu, dataset$mu_tsd, dataset$sigma_tmu, dataset$sigma_tsd
            )
        }
    }
    if (transform) {
        which_transf <- which(grepl("log", rnames))
        if (any(grepl("spp_", rnames))) {
            which_transf <- c(which_transf, which(grepl("spp_", rnames)))
        }
        if (any(grepl("sigma", rnames))) {
            which_transf <- c(which_transf, which(grepl("sigma", rnames)))
        }
        precis_result$mean[which_transf] <- exp(precis_result$mean[which_transf])
        precis_result$sd[which_transf] <- exp(precis_result$sd[which_transf])
        precis_result$`5.5%`[which_transf] <- exp(precis_result$`5.5%`[which_transf])
        precis_result$`94.5%`[which_transf] <- exp(precis_result$`94.5%`[which_transf])
    }
    precis_result$delta <- precis_result$expected - precis_result$mean
    if (change_names) {
        precis_result$what <- gsub("log_tmu", "tmu", precis_result$what)
        precis_result$what <- gsub("log_tsd", "tsd", precis_result$what)
    }
    return(precis_result)
}

sim_dataset_phylo <- function(
    N = 500,                                    # Number of surveys
    n_sp = 5,                                   # Number of species
    p = rbeta(1, 2, 2),                         # Detection probability
    global_optimum = 1.5,
    rho_val = 1,
    max_covariance = 3,
    equal_n_pres = FALSE,
    mu_tmu = 20,                                # Global mean for tmu
    sigma_tmu = 3,                              # Global sd for tmu
    mu_tsd = 5,                                 # Global mean for tsd
    sigma_tsd = 1,                              # Global sd for tsd 
    tomax = rbeta(n_sp, 5, 1),                  # Max occupancy prob per species - you can also pass a fixed value
    site_min = 5,                               # Minimum of site
    site_max = 30,                              # Maximum of site
    n_presence = NULL,                          # Number of presences
    n_absence = NULL,                           # Number of absences
    perc_sst_inc = 0,                           # Percentage of increase in SST SD
    make_plot = TRUE                            # Plot objects
) {
    require(ape)
    require(rethinking)

    log_inv_logit <- function(x) log(inv_logit(x))

    D <- 3 # dimension of random effects

    # simulate tree and compute distance matrix
    a_tree <- rtree(n_sp)
    M <- cophenetic(a_tree) # pairwise cophenetic distance matrix
    M <- M / max(M) # rescale dists to be 0/1
    Omega <- rlkjcorr(1, D, 1) # correlation structure of outcome vector

    if (make_plot) {
        par(mfrow = c(1, 2))
        plot(a_tree)
        image(M)
    }

    # build covariance kernel across spp
    K <- matrix(NA, n_sp, n_sp)
    rho <- rho_val
    # higher rho = more local correlation -> less overall independent spread for distant species
    # lower rho = longer-range correlation -> species remain similar over larger distances
    # If very high, than phylogeny no longer has an effect

    # L1 distance kernel
    max_cov <- max_covariance # In the phylogenetic signal => high closely related similar
    # Plot M = horizonal
    for (i in 1:n_sp) {
        for (j in 1:n_sp) {
            K[i, j] <- max_cov * exp(-rho * M[i, j]) # covariance decays with distance
        }
    }
    diag(K) <- max_cov + 0.001
    LK <- t(chol(K))
    LO <- t(chol(Omega))

    # we want f = LK * eta * LO'
    k <- c(global_optimum, 1, 2) # intercept
    # k[1] is the global mean for optimum temperature (log scale)
    f <- matrix(NA, n_sp, D)
    eta <- rmvnorm2(n_sp, sigma = rep(1, D), Rho = diag(D))
    f <- LK %*% eta %*% t(LO) # generates a n_sp x D matrix of latent values
    for (i in 1:n_sp) {
        f[i, ] <- f[i, ] + k # sum intercept
        if (exp(f[i, 1]) > 35) f[i, 1] <- log(35)
    }

    n_each <- N / n_sp
    sid <- rep(seq_len(n_sp), each = n_each)
    sst <- rep(NA, N)
    y <- rep(NA, N)

    # plot f scales
    if (make_plot) {
        par(mfrow = c(1, 2))
        plot(exp(f[sid, 1]), exp(f[sid, 2]))
        plot(exp(f[sid, 1]), inv_logit(f[sid, 3]))
    }

    if (equal_n_pres) {
        pts_p_sp <- N/n_sp
        for (idsp in seq_len(n_sp)) {
            presences <- floor(pts_p_sp / 2)
            total_p <- 0
            total_a <- 0
            max_try <- 10000
            tried <- 0
            sst_sp <- rep(NA, max_try)
            y_sp <- rep(NA, max_try)
            idx <- 1
            samp_pts <- function(sid) {
                sst <- runif(1, site_min, site_max) # temperature at site
                log_q <- (-0.5) * ((sst - exp(f[sid, 1])) / exp(f[sid, 2]))^2 + log_inv_logit(f[sid, 3])
                Y <- rbern(1, exp(log_q))
                y <- rbern(1, Y * p) # With detection probability
                list(sst = sst, y = y)
            }
            while((total_p < presences || total_a < presences) & tried < max_try) {
                sr <- samp_pts(idsp)
                sst_sp[idx] <- sr$sst
                y_sp[idx] <- sr$y
                total_p <- sum(y_sp, na.rm = TRUE)
                total_a <- sum(y_sp == 0, na.rm = TRUE)
                idx <- idx + 1
                tried <- tried + 1
            }
            if (total_p < presences || total_a < presences) {
                warning("Not possible to get equal numbers for species ", idsp)
                sst[sid == idsp] <- sst_sp[seq_len(presences * 2)]
                y[sid == idsp] <- y_sp[seq_len(presences * 2)]
            } else {
                samp_p <- sample(which(y_sp == 1), presences)
                samp_a <- sample(which(y_sp == 0), presences)
                sst[sid == idsp] <- sst_sp[c(samp_p, samp_a)]
                y[sid == idsp] <- y_sp[c(samp_p, samp_a)]
            }
        }
    } else {
        for (i in seq_len(N)) {
            sst[i] <- runif(1, site_min, site_max) # temperature at site
            # log prob occupy site
            # log_q <- [log of a normal density => center = exp(f[sid[i], 1])), spread =  exp(f[sid[i], 2])] + log_inv_logit(f3)
            log_q <- (-0.5) * ((sst[i] - exp(f[sid[i], 1])) / exp(f[sid[i], 2]))^2 + log_inv_logit(f[sid[i], 3])
            # sim true occupancy
            Y <- rbern(1, exp(log_q))
            # sim observation
            y[i] <- rbern(1, Y * p) # With detection probability
        }
    }

    
    # species_params <- data.frame(
    # species = 1:nrow(f),
    # opt_temp = exp(f[, 1]),
    # temp_width = exp(f[, 2]),
    # occupancy_prob = plogis(f[, 3])
    # )
    # print(species_params)

    dat_res <- list(
        data = list(
            N = N,
            D = D,
            N_spp = n_sp,
            sid = sid,
            sst = sst,
            y = y,
            M = M
        ),
        raw = list(
            f = f
        )
    )

    return(dat_res)
}
