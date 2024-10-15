
dco <- function(community_matrix, mcmc_control) {
    # check control
    mcmc_control <- process_mcmc_control(control_in = mcmc_control)

    if(ncol(community_matrix)==1) {
        stop("requirement for dco analysis not met: ncol(community_matrix)>1")
    }
    if(ncol(community_matrix)==2) {
        model <- stanmodels$dm_n2
    }
    if(ncol(community_matrix)>2) {
        model <- stanmodels$dm
    }
    
    # fit
    message("1/2 fit...")
    f <- sampling(object = model,
                  data = list(K = nrow(community_matrix), 
                              N = ncol(community_matrix), 
                              y = t(community_matrix),
                              x = c(-1, +1)),
                  chains = mcmc_control$mcmc_chains, 
                  cores = mcmc_control$mcmc_cores, 
                  iter = mcmc_control$mcmc_iter, 
                  warmup = mcmc_control$mcmc_warmup,
                  control = list(adapt_delta = mcmc_control$adapt_delta,
                                 max_treedepth = mcmc_control$max_treedepth),
                  algorithm = mcmc_control$mcmc_algorithm,
                  include = TRUE,
                  pars = c("alpha", "beta", "kappa", "p", "y_hat"))
    # summaries
    message("2/2 posterior summary...")
    s <- get_posterior_summaries(cm = community_matrix, f = f)
    
    return(list(fit = f, 
                posterior_summary = s, 
                community_matrix = community_matrix, 
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
    
    get_sample_com_par <- function(f, samples, par) {
        s <- data.frame(summary(f, par = par)$summary)
        
        # maintain original index order
        s$i <- 1:nrow(s)
        browser()
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
            s <- merge(x = s, y = meta, by.x = "s", by.y = "s", all.x = T)
            s <- s[order(s$i, decreasing = F),]
        }
        
        return(s)
    }
    
    get_com_par <- function(f, par) {
        s <- data.frame(summary(f, par = par)$summary)
        
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        s$community <- as.numeric(m)
        return(s)
    }
    
    get_sample_par <- function(f, par) {
        s <- data.frame(summary(f, par = par)$summary)
        
        m <- rownames(s)
        m <- gsub(pattern = paste0(par,"\\[|\\]"), replacement = '', x = m)
        
        s$s <- as.numeric(m)
        meta <- data.frame(s = 1:length(samples), sample = samples)
        s <- merge(x = s, y = meta, by.x = "s", by.y = "s", all.x = T)
        
        return(s)
    }
    
    get_global_par <- function(f, par) {
        s <- data.frame(summary(f, par = par)$summary)
        s$par <- rownames(s)
        return(s)
    }

    samples <- colnames(cm)
    
    o <- list(beta = get_sample_com_par(f = f, samples = samples, par = "beta"),
              alpha = get_com_par(f = f, par = "alpha"),
              p = get_sample_com_par(f = f, samples = samples, par = "p"),
              y_hat=get_sample_com_par(f = f, samples = samples, par = "y_hat"),
              kappa = get_global_par(f = f, par = "kappa"))
    
    o$y_hat$y_obs <- c(cm)
    
    return(o)
}

