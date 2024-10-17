
dco <- function(community_occupancy_matrix, mcmc_control) {
    # check control
    mcmc_control <- process_mcmc_control(control_in = mcmc_control)

    n <- ncol(community_occupancy_matrix)
    if(n==1) {
        stop("ncol(community_occupancy_matrix) must be >1")
    }
    if(n==2) {
        model <- stanmodels$dm_n2
        pars <- c("alpha", "beta", "kappa", "p", "y_hat", "log_lik")
    }
    if(n>2) {
        model <- stanmodels$dm
        pars <- c("alpha", "beta", "delta", "kappa", "p", "y_hat", "log_lik")
    }
    
    # fit
    message("1/2 fit...")
    f <- sampling(object = model,
                  data = list(K = nrow(community_occupancy_matrix), 
                              N = ncol(community_occupancy_matrix), 
                              y = t(community_occupancy_matrix),
                              x = c(-1, +1)),
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
    message("2/2 posterior summary...")
    s <- get_posterior_summaries(cm = community_occupancy_matrix, f = f)
    
    return(list(fit = f, 
                posterior_summary = s, 
                community_occupancy_matrix = community_occupancy_matrix, 
                mcmc_control = mcmc_control))
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


get_posterior_summaries <- function(cm, f) {
    
    post_sample_com <- function(f, samples, par) {
        
        s <- data.frame(summary(f, par = par)$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        # maintain original index order
        s$i <- 1:nrow(s)
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        
        if(length(samples)==2 & par == "beta") {
            s$community <- as.numeric(m)
        } 
        else {
            m <- do.call(rbind, strsplit(x = m, split = "\\,"))
            s$s <- as.numeric(m[,1])
            s$community <- as.numeric(m[,2])
            meta <- data.frame(s = 1:length(samples), sample = samples)
            s <- merge(x = s, y = meta, by.x = "s", by.y = "s", all.x = TRUE)
            s <- s[order(s$i, decreasing = FALSE),]
        }
        
        return(s)
    }
    
    post_com <- function(f, par) {
        
        s <- data.frame(summary(f, par = par)$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        s$community <- as.numeric(m)
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
    
    post_global <- function(f, par) {
        
        s <- data.frame(summary(f, par = par)$summary)
        s <- s[, c("mean", "X50.", "X2.5.", "X97.5.", "n_eff", "Rhat")]
        colnames(s) <- c("mean", "median", "L95", "H95", "n_eff", "Rhat")
        
        s$par <- rownames(s)
        return(s)
    }
    
    post_delta <- function(f, samples, par) {
        
        if(length(samples)==2) {
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
        
        s <- merge(x = s, y = meta, by.x = "k", by.y = "k", all.x = T)
        s <- s[order(s$i, decreasing = F),]
        s$i <- NULL
        
        # get reverse effects
        sr <- s
        sr$mean <- sr$mean*-1
        sr$median <- sr$median*-1
        sr$L95 <- sr$L95*-1
        sr$H95 <- sr$H95*-1
        x <- sr$sample_1
        y <- sr$sample_2
        sr$sample_2 <- x
        sr$sample_1 <- y
        
        s <- rbind(s, sr)
        s$contrast <- paste0(s$sample_1, '-', s$sample_2)
        
        return(s)
    }

    samples <- colnames(cm)
    
    o <- list(beta = post_sample_com(f = f, samples = samples, par = "beta"),
              alpha = post_com(f = f, par = "alpha"),
              p = post_sample_com(f = f, samples = samples, par = "p"),
              y_hat = post_sample_com(f = f, samples = samples, par = "y_hat"),
              kappa = post_global(f = f, par = "kappa"),
              delta = post_delta(f = f, samples = samples, par = "delta"))
    
    o$y_hat$y_obs <- c(cm)

    return(o)
}

