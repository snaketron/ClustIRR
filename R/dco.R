
dco <- function(community_occupancy_matrix, 
                mcmc_control, 
                compute_delta = TRUE,
                groups = NA) {
    # check control
    mcmc_control <- process_mcmc_control(control_in = mcmc_control)
    
    n <- ncol(community_occupancy_matrix)
    if(n==1) {
        stop("ncol(community_occupancy_matrix) must be >1")
    }
    
    # check compute_delta
    if(missing(compute_delta)) {
        compute_delta <- TRUE
    }
    if(length(compute_delta)!=1) {
        stop("compute_delta must be logical")
    }
    if(is.logical(compute_delta)==FALSE) {
        stop("compute_delta must be logical")
    }
    compute_delta <- ifelse(test = compute_delta == TRUE, yes = 1, no = 0)
    
    
    # check groups
    if(missing(groups)) {
        groups <- NA
    }
    has_groups <- all(is.na(groups)==FALSE)
    if(has_groups) {
        if(length(groups)!=n) {
            stop("length(groups) != ncol(community_occupancy_matrix)")
        }
        if(all(is.numeric(groups))==FALSE) {
            stop("groups must integers")
        }
        if(all(1:max(groups) %in% groups)==FALSE) {
            stop("missing group in 1:max(groups)")
        }
    }
    
    # select normal or hierarchical model
    if(has_groups) {
        model <- stanmodels$dmh
        pars <- c("alpha", "beta", "beta_mu", "beta_sigma", 
                  "kappa", "y_hat", "log_lik")
        if(compute_delta == 1) {
            pars <- c("alpha", "beta", "beta_mu", "beta_sigma", 
                      "delta", "epsilon", "kappa", "y_hat", "p", "log_lik")
        }
    } 
    else {
        model <- stanmodels$dm
        pars <- c("alpha", "beta", "kappa", "y_hat", "log_lik")
        if(compute_delta == 1) {
            pars <- c("alpha", "beta", "delta", "epsilon", "kappa", 
                      "y_hat", "p", "log_lik")
        }
    }
    
    
    # fit
    message("[1/2] fit...")
    f <- sampling(object = model,
                  data = list(K = nrow(community_occupancy_matrix), 
                              N = n, 
                              y = t(community_occupancy_matrix),
                              G = groups,
                              compute_delta = compute_delta),
                  chains = mcmc_control$mcmc_chains, 
                  cores = mcmc_control$mcmc_cores, 
                  iter = mcmc_control$mcmc_iter, 
                  warmup = mcmc_control$mcmc_warmup,
                  control = list(adapt_delta = mcmc_control$adapt_delta,
                                 max_treedepth = mcmc_control$max_treedepth),
                  algorithm = mcmc_control$mcmc_algorithm,
                  include = TRUE,
                  pars = pars)
    
    # summaries
    message("[2/2] posterior summary...")
    s <- get_posterior_summaries(cm = community_occupancy_matrix, 
                                 f = f, 
                                 groups = groups,
                                 has_groups = has_groups,
                                 compute_delta = compute_delta)
    
    return(list(fit = f, 
                posterior_summary = s, 
                community_occupancy_matrix = community_occupancy_matrix, 
                mcmc_control = mcmc_control,
                compute_delta = compute_delta,
                groups = groups))
}


process_mcmc_control <- function(control_in) {
    
    check_mcmc_iter <- function(mcmc_iter,
                                mcmc_warmup) {
        
        if(length(mcmc_iter) != 1) {
            stop("mcmc_iter must be a positive integer")
        }
        if(length(mcmc_warmup) != 1) {
            stop("mcmc_warmup must be a positive integer")
        }
        
        if(is.numeric(mcmc_iter) == FALSE) {
            stop("mcmc_iter must be a positive integer")
        }
        if(is.numeric(mcmc_warmup) == FALSE) {
            stop("mcmc_warmup must be a positive integer")
        }
        
        if(is.finite(mcmc_iter) == FALSE) {
            stop("mcmc_iter must be a positive integer")
        }
        if(is.finite(mcmc_warmup) == FALSE) {
            stop("mcmc_warmup must be a positive integer")
        }
        
        if((abs(mcmc_iter - round(mcmc_iter)) 
            < .Machine$double.eps^0.5) == FALSE) {
            stop("mcmc_iter must be a positive integer")
        }
        if((abs(mcmc_warmup - round(mcmc_warmup)) 
            < .Machine$double.eps^0.5) == FALSE) {
            stop("mcmc_warmup must be a positive integer")
        }
        
        
        if(mcmc_iter <= mcmc_warmup) {
            stop("mcmc_iter must be larger than mcmc_warmup")
        }
    }
    
    check_mcmc_chains <- function(mcmc_chains) {
        if(length(mcmc_chains) != 1) {
            stop("mcmc_chains must be a positive integer")
        }
        
        if(!is.numeric(mcmc_chains)) {
            stop("mcmc_chains must be a positive integer")
        }
        
        if(is.finite(x = mcmc_chains) == FALSE) {
            stop("mcmc_chains must be a positive integer")
        }
        
        if(as.integer(x = mcmc_chains) <= 0) {
            stop("mcmc_chains must be a positive integer")
        }
    }
    
    check_mcmc_cores <- function(mcmc_cores) {
        if(length(mcmc_cores) != 1) {
            stop("mcmc_cores must be a positive integer")
        }
        
        if(is.numeric(mcmc_cores) == FALSE) {
            stop("mcmc_cores must be a positive integer")
        }
        
        if(is.finite(x = mcmc_cores) == FALSE) {
            stop("mcmc_cores must be a positive integer")
        }
        
        if(as.integer(x = mcmc_cores) <= 0) {
            stop("mcmc_cores must be a positive integer")
        }
    }
    
    check_adapt_delta <- function(adapt_delta) {
        if(length(adapt_delta) != 1) {
            stop("adapt_delta must be in the range between 0 and 1")
        }
        
        if(is.numeric(adapt_delta) == FALSE) {
            stop("adapt_delta must be in the range between 0 and 1")
        }
        
        if(is.finite(x = adapt_delta) == FALSE) {
            stop("adapt_delta must be in the range between 0 and 1")
        }
        
        if(adapt_delta<0 | adapt_delta>1) {
            stop("adapt_delta must be in the range between 0 and 1")
        }
    }
    
    check_max_treedepth <- function(max_treedepth) {
        if(length(max_treedepth) != 1) {
            stop("max_treedepth must be a positive integer")
        }
        
        if(is.numeric(max_treedepth) == FALSE) {
            stop("max_treedepth must be a positive integer")
        }
        
        if(is.finite(x = max_treedepth) == FALSE) {
            stop("max_treedepth must be a positive integer")
        }
        
        if((abs(max_treedepth - round(max_treedepth)) 
            < .Machine$double.eps^0.5) == FALSE) {
            stop("max_treedepth must be a positive integer")
        }
        
        if(as.integer(x = max_treedepth) <= 0) {
            stop("max_treedepth must be a positive integer")
        }
    }
    
    check_mcmc_algorithm <- function(mcmc_algorithm) {
        if(length(mcmc_algorithm) != 1) {
            stop("mcmc_algorithm must be set to NUTS")
        }
        
        if(is.character(mcmc_algorithm) == FALSE) {
            stop("mcmc_algorithm must be set to NUTS")
        }
        
        if(mcmc_algorithm!="NUTS") {
            stop("mcmc_algorithm must be set to NUTS")
        }
    }
    
    
    control <- list(mcmc_warmup = 750,
                    mcmc_iter = 1500,
                    mcmc_chains = 4,
                    mcmc_cores = 1,
                    mcmc_algorithm = "NUTS",
                    adapt_delta = 0.95,
                    max_treedepth = 12)
    
    # if missing control_in -> use default values
    if(missing(control_in) || is.null(control_in)) {
        return(control)
    }
    if(is.list(control_in) == FALSE) {
        stop("control must be a list")
    }
    if(all(names(control_in) %in% names(control)) == FALSE) {
        stop("unrecognized elements found in control")
    }
    
    ns <- names(control_in)
    for (i in seq_len(length(control_in))) {
        control[[ns[i]]] <- control_in[[ns[i]]]
    }
    
    
    check_mcmc_iter(mcmc_iter = control$mcmc_iter,
                    mcmc_warmup = control$mcmc_warmup)
    check_mcmc_chains(mcmc_chains = control$mcmc_chains)
    check_mcmc_cores(mcmc_cores = control$mcmc_cores)
    check_adapt_delta(adapt_delta = control$adapt_delta)
    check_max_treedepth(max_treedepth = control$max_treedepth)
    check_mcmc_algorithm(mcmc_algorithm = control$mcmc_algorithm)
    return(control)
}


get_posterior_summaries <- function(cm, 
                                    f, 
                                    compute_delta, 
                                    has_groups, 
                                    groups) {
    
    post_yhat <- function(f, samples) {
        
        s <- data.frame(summary(f, par = "y_hat")$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")

        # maintain original index order
        s$i <- 1:nrow(s)
        m <- rownames(s)
        m <- gsub(pattern = "y_hat\\[|\\]", replacement = '', x = m)
        
        m <- do.call(rbind, strsplit(x = m, split = "\\,"))
        s$s <- as.numeric(m[,1])
        s$community <- as.numeric(m[,2])
        meta <- data.frame(s = 1:length(samples), sample = samples)
        s <- merge(x = s, y = meta, by.x = "s", by.y = "s", all.x = TRUE)
        s <- s[order(s$i, decreasing = FALSE),]
        
        return(s)
    }
    
    post_alpha <- function(f) {
        
        s <- data.frame(summary(f, par = "alpha")$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        m <- rownames(s)
        m <- gsub(pattern = paste0("alpha\\[|\\]"), replacement = '', x = m)
        s$community <- as.numeric(m)
        return(s)
    }
    
    post_beta <- function(f, samples) {
        
        if(length(samples) == 2) {
            s <- data.frame(summary(f, par = "beta")$summary)
            s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
            colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
            
            s$i <- 1:nrow(s)
            m <- rownames(s)
            m <- gsub(pattern = paste0("beta\\[|\\]"), replacement = '', x = m)
            m <- do.call(rbind, strsplit(x = m, split = "\\,"))
            s$s <- as.numeric(m[,1])
            s$community <- as.numeric(m[,2])
            s <- s[order(s$i, decreasing = FALSE),]
        } 
        else {
            s <- data.frame(summary(f, par = "beta")$summary)
            s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
            colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
            
            # maintain original index order
            s$i <- 1:nrow(s)
            m <- rownames(s)
            m <- gsub(pattern = "beta\\[|\\]", replacement = '', x = m)
            
            m <- do.call(rbind, strsplit(x = m, split = "\\,"))
            s$s <- as.numeric(m[,1])
            s$community <- as.numeric(m[,2])
            meta <- data.frame(s = 1:length(samples), sample = samples)
            s <- merge(x = s, y = meta, by.x = "s", by.y = "s", all.x = TRUE)
            s <- s[order(s$i, decreasing = FALSE),]
        }
        
        return(s)
    }
    
    post_sample <- function(f, par) {
        
        s <- data.frame(summary(f, par = par)$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        
        s$s <- as.numeric(m)
        meta <- data.frame(s = 1:length(samples), sample = samples)
        s <- merge(x = s, y = meta, by.x = "s", by.y = "s", all.x = TRUE)
        
        return(s)
    }
    
    post_kappa <- function(f) {
        s <- data.frame(summary(f, par = "kappa")$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        s$par <- rownames(s)
        return(s)
    }
    
    post_betamu <- function(f, has_groups) {
        if(has_groups==FALSE) {
            return(NA)
        }
        
        s <- data.frame(summary(f, par = "beta_mu")$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        # maintain original index order
        s$i <- 1:nrow(s)
        m <- rownames(s)
        m <- gsub(pattern = "beta\\_mu\\[|\\]", replacement = '', x = m)
        
        m <- do.call(rbind, strsplit(x = m, split = "\\,"))
        s$g <- as.numeric(m[,1])
        s$community <- as.numeric(m[,2])
        
        return(s)
    }
    
    post_betasigma <- function(f, has_groups) {
        if(has_groups==FALSE) {
            return(NA)
        }
        s <- data.frame(summary(f, par = "beta_sigma")$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        m <- rownames(s)
        m <- gsub(pattern = "beta\\_sigma\\[|\\]", replacement = '', x = m)
        s$group <- as.numeric(m)
        return(s)
    }
    
    post_delta <- function(f, samples, compute_delta, par) {
        if(compute_delta==0) {
            return(NA)
        }
        
        get_meta <- function(samples) {
            m <- c()
            k <- 1
            for(i in 1:(length(samples)-1)) {
                for(j in (i+1):length(samples)) {
                    m <- rbind(m, data.frame(k = k, 
                                             sample_1 = samples[i],
                                             sample_2 = samples[j]))
                    k <- k + 1
                }   
            }
            m$contrast <- paste0(m$sample_1, '-', m$sample_2)
            return(m)
        }
        
        meta <- get_meta(samples = samples)
        
        # pmax
        p <- as_draws(f)
        p <- subset_draws(p, variable = par)
        p <- summarise_draws(p, pmax =~2*max(mean(.>0), mean(.<=0))-1)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', 
                  x = p$variable)
        m <- do.call(rbind, strsplit(x = m, split = "\\,"))
        p$k <- as.numeric(m[,1])
        p$community <- as.numeric(m[,2])
        
        # main summary
        s <- data.frame(summary(f, par = par)$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        # maintain original index order
        s$i <- 1:nrow(s)
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        m <- do.call(rbind, strsplit(x = m, split = "\\,"))
        s$k <- as.numeric(m[,1])
        s$community <- as.numeric(m[,2])
        
        
        s <- merge(x = s, y = p[,c("k", "community", "pmax")], 
                   by = c("k", "community"))
        s <- merge(x = s, y = meta, by.x = "k", by.y = "k", all.x = TRUE)
        s <- s[order(s$i, decreasing = FALSE),]
        s$i <- NULL
        
        
        # get reverse effects
        sr <- s
        sr$mean <- sr$mean*-1
        sr$median <- sr$median*-1
        x <- sr$L95*-1
        y <- sr$H95*-1
        sr$L95 <- y
        sr$H95 <- x
        x <- sr$sample_1
        y <- sr$sample_2
        sr$sample_2 <- x
        sr$sample_1 <- y
        
        s <- rbind(s, sr)
        s$contrast <- paste0(s$sample_1, '-', s$sample_2)
        
        return(s)
    }
    
    post_deltamu <- function(f, groups, compute_delta, par) {
        if(compute_delta==0) {
            return(NA)
        }
        
        get_meta <- function(groups) {
            m <- c()
            k <- 1
            for(i in 1:(length(groups)-1)) {
                for(j in (i+1):length(groups)) {
                    m <- rbind(m, data.frame(k = k, 
                                             group_1 = groups[i],
                                             group_2 = groups[j]))
                    k <- k + 1
                }   
            }
            m$contrast <- paste0(m$group_1, '-', m$group_2)
            return(m)
        }
        
        meta <- get_meta(groups = groups)
        
        # pmax
        p <- as_draws(f)
        p <- subset_draws(p, variable = par)
        p <- summarise_draws(p, pmax =~2*max(mean(.>0), mean(.<=0))-1)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = p$variable)
        m <- do.call(rbind, strsplit(x = m, split = "\\,"))
        p$k <- as.numeric(m[,1])
        p$community <- as.numeric(m[,2])
        
        # main summary
        s <- data.frame(summary(f, par = par)$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        # maintain original index order
        s$i <- 1:nrow(s)
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        m <- do.call(rbind, strsplit(x = m, split = "\\,"))
        s$k <- as.numeric(m[,1])
        s$community <- as.numeric(m[,2])
        
        
        s <- merge(x = s, y = p[,c("k", "community", "pmax")], 
                   by = c("k", "community"))
        s <- merge(x = s, y = meta, by.x = "k", by.y = "k", all.x = TRUE)
        s <- s[order(s$i, decreasing = FALSE),]
        s$i <- NULL
        
        
        # get reverse effects
        sr <- s
        sr$mean <- sr$mean*-1
        sr$median <- sr$median*-1
        x <- sr$L95*-1
        y <- sr$H95*-1
        sr$L95 <- y
        sr$H95 <- x
        x <- sr$group_1
        y <- sr$group_2
        sr$group_2 <- x
        sr$group_1 <- y
        
        s <- rbind(s, sr)
        s$contrast <- paste0(s$group_1, '-', s$group_2)
        
        return(s)
    }
    
    samples <- colnames(cm)
    groups <- unique(groups)
    
    o <- list(beta = post_beta(f = f, samples = samples),
              beta_mu = post_betamu(f = f, has_groups = has_groups),
              beta_sigma = post_betasigma(f = f, has_groups = has_groups),
              alpha = post_alpha(f = f),
              y_hat = post_yhat(f = f, samples = samples),
              kappa = post_kappa(f = f))
    
    # compute delta
    if(has_groups==FALSE) {
        delta <- post_delta(f = f, samples = samples, 
                            compute_delta = compute_delta, par = "delta")
        epsilon <- post_delta(f = f, samples = samples,
                              compute_delta = compute_delta, par = "epsilon")
    } else {
        delta <- post_deltamu(f = f, groups = groups, 
                              compute_delta = compute_delta, par = "delta")
        epsilon <- post_deltamu(f = f, groups = groups,
                                compute_delta = compute_delta, par = "epsilon")
    }
    o[["delta"]] <- delta
    o[["epsilon"]] <- epsilon
    o$y_hat$y_obs <- c(cm)
    
    return(o)
}
