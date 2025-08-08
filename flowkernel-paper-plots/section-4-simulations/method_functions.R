# method_functions.R

kernel_em_method <- new_method(name = "kernel_em",
                               label = "Kernel EM Fit",
                               method = function(model, draw) {
                                 #crossprod <- crossprod
                                 # Extract data from the draw
                                 dat <- draw
                                 K <- 2
                                 hmu <- model$params$hmu
                                 hSigma <- model$params$hSigma
                                 hpi <- model$params$hpi
                                 
                                 # Initialize and fit the model
                                 initial_fit <- init_const(dat$y, K = K)
                                 fit <- kernel_em_dates(dat$y, K = K, hmu = hmu, hSigma = hSigma, hpi = hpi,
                                                        dates = dat$dates, biomass = dat$biomass,
                                                        initial_fit = initial_fit)
                                 
                                 # Return the fit and true labels
                                 list(fit = fit, true_z = dat$z)
                               })

kernel_em_CV_method <- new_method(name = "kernel_em_cv",
                                  label = "Kernel EM Fit",
                                  method = function(model, draw) {
                  
                                    # Extract data from the draw
                                    dat <- draw
                                    K <- 2
                                    hSigma <- 100
                                    #disappearance_duration <- as.character(model$params$disappearance_duration)
                                    
                                    #cv_results_par <- cv_log_grid_search(data = dat$y, K = K, hSigma = hSigma,dates = dat$dates, biomass = dat$biomass,leave_out_every = 5, grid_size = 7)
                                    cv_results_par <- cv_log_grid_search_parallel(data = dat$y, K = K, hSigma = hSigma,
                                                                                  dates = dat$dates, biomass = dat$biomass,
                                                                                  leave_out_every = 5, grid_size = 7)
                                    
                                    if (is.null(cv_results_par$best_hmu) || is.null(cv_results_par$best_hpi)) {
                                      stop("Error: Cross-validation did not return a valid best_hmu or best_hpi.")
                                    }
                                    hmu <- cv_results_par$best_hmu
                                    hpi <- cv_results_par$best_hpi

                                    
                                    # Initialize and fit the model
                                    initial_fit <- init_const(dat$y, K = K)
                                    fit <- kernel_em_dates(dat$y, K = K, hmu = hmu, hSigma = hSigma, hpi = hpi,
                                                           dates = dat$dates, biomass = dat$biomass,
                                                           initial_fit = initial_fit)
                                    
                                    # Return the fit and true labels
                                    list(fit = fit, true_z = dat$z)
                                  })



const_method <- new_method(name = "constant_fit",
                           label = "Constant Fit",
                           method = function(model, draw) {
                             #crossprod <- crossprod
                             # Extract data from the draw
                             dat <- draw
                             K <- 2
                             # Initialize and fit the model
                             fit <- const_fit(dat$y, K = K)
                             # Return the fit and true labels
                             list(fit = fit, true_z = dat$z)
                           })


library(mclust)
library(RcppHungarian)
hungarian_method <- new_method(name = "hungarian_fit",
                               label = "Hungarian Fit",
                               method = function(model, draw) {
                                 #crossprod <- crossprod
                                 # Extract data from the draw
                                 dat <- draw
                                 K <- 2
                                 # fit the model
                                 fit <- hungarian_fit(dat$y, K = K)
                                 
                                 # Return the fit and true labels
                                 list(fit = fit, true_z = dat$z)
                               })

hungarian_fit <- function(y, K) {
  num_times <- length(y)
  d <- ncol(y[[1]])
  mu <- array(NA, c(num_times, K, d))
  Sigma <- array(NA, c(num_times, K, d, d))
  pi <- matrix(NA, num_times, K)
  fit_save <- NULL  # To save the last successful fit
  
  if (d == 1) {
    for (tt in seq_len(num_times)) {
      print(tt)
      init_fit <- mclust::Mclust(y[[tt]], G = K, modelNames = "V")
      if (is.null(init_fit) || init_fit$G != K) {
        # If we have a valid 'fit_save' from a previous iteration that had G=2, use it.
        if (!is.null(fit_save) && !is.null(fit_save$parameters) && fit_save$G == 2) {
          init_fit <- fit_save
        } else {
          # We have no valid 2-cluster fallback. 
          # Forcibly invent a 2-cluster assignment.
          mu[tt, , 1] <- rep(mean(y[[tt]]), K)
          Sigma[tt, , 1, 1] <- rep(var(y[[tt]]), K)
          pi[tt, ] <- c(0.5, 0.5)
          # Then move on to the next time step, or break, or whatever is appropriate.
          next
        }
      }
      fit_save <- init_fit
      mu[tt, , 1] <- init_fit$parameters$mean
      Sigma[tt, , 1, 1] <- init_fit$parameters$variance$sigmasq
      pi[tt, ] <- init_fit$parameters$pro
    }
  } else if (d > 1) {  # Multivariate initialization
    for (tt in seq_len(num_times)) {
      print(tt)
      init_fit <- mclust::Mclust(y[[tt]], G = K, modelNames = "VVV")
      if (is.null(init_fit)) {
        init_fit <- fit_save
        print(tt * 1000)
      }
      fit_save <- init_fit
      mu[tt, , ] <- t(init_fit$parameters$mean)
      pi[tt, ] <- init_fit$parameters$pro
      Sigma[tt, , , ] <- aperm(init_fit$parameters$variance$sigma, c(3, 1, 2))
    }
    
    # Use Hungarian Algorithm to ensure consistent labeling of clusters over time
    for (t in 2:num_times) {
      permu_vec <- permute_cluster(mu[t - 1, , ], mu[t, , ])
      mu[t, , ] <- mu[t, , ][permu_vec, ]
      pi[t, ] <- pi[t, ][permu_vec]
      Sigma[t, , , ] <- Sigma[t, , , ][permu_vec, , ]
    }
  }
  resp <- calculate_responsibilities(y, mu, Sigma, pi)
  zest <- resp %>% purrr::map(~ max.col(.x))
  list(mu = mu, Sigma = Sigma, pi = pi, resp = resp, zest = zest)
}

cost_dist <- function(mu1, mu2){
  K <- nrow(mu1)
  C = matrix(NA, nrow = K, ncol = K)
  for (i in seq(K)){
    for (j in seq(K)){
      C[i,j] = dist(rbind(mu1[i,], mu2[j,]))
    }
  }
  return(C)
}

permute_cluster = function(mu1,mu2) {#mu1 has its clusters fixed, mu2 will be rearranged
  cost_mat = cost_dist(mu1,mu2)
  permu_vec = HungarianSolver(cost_mat)$pairs[,2]
  return(permu_vec)
}

const_fit <- function (y, K) {
  num_times <- length(y)
  d <- ncol(y[[1]])
  mu <- array(NA, c(num_times, K, d))
  Sigma <- array(NA, c(num_times, K, d, d))
  pi <- matrix(NA, num_times, K)
  
  # Repeatedly call Mclust until it gives a non-NULL fit
  init_fit <- NULL
  count <- 1
  while (is.null(init_fit)) {
    # Combine data from all times into a single matrix
    y_combined <- do.call(rbind, y)
    
    if (count == 1) {
      y_sampled <- y_combined
    } else {
      # Sample 90% of the data
      sample_indices <- sample(1:nrow(y_combined), size = floor(0.9 * nrow(y_combined)))
      y_sampled <- y_combined[sample_indices, , drop = FALSE]
    }
    count <- count + 1
    
    # Try fitting Mclust on the sampled data
    if (d == 1) {
      init_fit <- mclust::Mclust(y_sampled, G = K, modelNames = "V")
      if (!is.null(init_fit)) {
        for (tt in seq(num_times)) {
          mu[tt, , 1] <- init_fit$parameters$mean
          Sigma[tt, , 1, 1] <- init_fit$parameters$variance$sigmasq
          pi[tt, ] <- init_fit$parameters$pro
        }
      }
    } else if (d > 1) {
      init_fit <- mclust::Mclust(y_sampled, G = K, modelNames = "VVV")
      if (!is.null(init_fit) && is.matrix(init_fit$parameters$mean)) {
        for (tt in seq(num_times)) {
          mu[tt, , ] <- t(init_fit$parameters$mean)
          pi[tt, ] <- init_fit$parameters$pro
          Sigma[tt, , , ] <- aperm(init_fit$parameters$variance$sigma, c(3, 1, 2))
        }
      }
    }
  }
  
  # Calculate responsibilities
  resp <- calculate_responsibilities(y, mu, Sigma, pi)
  zest <- resp %>% purrr::map(~ max.col(.x))
  
  list(mu = mu, Sigma = Sigma, pi = pi, resp = resp, zest = zest)
}

init_bayes <- function (y, K,  biomass = default_biomass(y), num_iter = 1,
                        lap_smooth_const = 0){
  num_times <- length(y)
  d <- ncol(y[[1]])
  resp <- list() # responsibilities gamma[[t]][i, k]
  log_resp <- list()
  resp_weighted <- list() 
  resp_sum <- list() 
  resp_sum_pi <- list()
  y_sum <- list()
  mat_sum <- list()
  yy <- list()
  mu <- array(NA, c(num_times, K, d))
  Sigma <- array(NA, c(num_times, K, d, d))
  pi <- matrix(NA, num_times, K)
  
  if (d == 1){
    init_fit <- mclust::Mclust(y[[1]], G = K, modelNames = "V")
    ii = 1
    while (is.null(init_fit) == TRUE){
      init_fit <- mclust::Mclust(y[[1 + ii]], G = K, modelNames = "V")
      ii = ii + 1
    }
    mu[1, , 1] <- init_fit$parameters$mean
    Sigma[1, , 1, 1] <- init_fit$parameters$variance$sigmasq
    pi [1, ] <- init_fit$parameters$pro
  } else if (d > 1){
    init_fit <- mclust::Mclust(y[[1]], G = K, modelNames = "VVV")
    ii = 1
    while (is.null(init_fit) == TRUE){
      init_fit <- mclust::Mclust(y[[1 + ii]], G = K, modelNames = "V")
      ii = ii + 1
    }
    mu[1, ,] <- t(init_fit$parameters$mean)
    pi[1, ] <- init_fit$parameters$pro
    Sigma[1, , , ] <- aperm(init_fit$parameters$variance$sigma, c(3,1,2))
  }
  phi <- matrix(NA, nrow(y[[1]]), K)
  log_phi <- matrix(NA, nrow(y[[1]]), K)
  if (d == 1) {
    for (k in seq(K)) {
      log_phi[, k] <- stats::dnorm(y[[1]],
                                   mean = mu[1, k, 1],
                                   sd = sqrt(Sigma[1, k, 1, 1]), log = TRUE)
    }
  } else if (d > 1) {
    for (k in seq(K)) {
      log_phi[, k] <- mvtnorm::dmvnorm(y[[1]],
                                       mean = mu[1, k, ],
                                       sigma = Sigma[1, k, , ], log = TRUE)
    }
  }
  log_temp = t(t(log_phi) + log(pi[1, ]))
  log_resp = log_temp - matrixStats::rowLogSumExps(log_temp)
  resp[[1]] = exp(log_resp)
  resp_weighted[[1]] = diag(biomass[[1]]) %*% resp[[1]]
  resp_sum[[1]] <- colSums(resp_weighted[[1]]) %>%
    unlist() %>%
    matrix(ncol = K, byrow = TRUE)
  pi_prev <- pi
  mu_prev <- mu
  Sigma_prev <- Sigma
  for (tt in 2:num_times){
    for (l in seq(num_iter)) {
      # E-step
      phi <- matrix(NA, nrow(y[[tt]]), K)
      if (l == 1){
        if (d == 1) {
          for (k in seq(K)) {
            phi[, k] <- stats::dnorm(y[[tt]],
                                     mean = mu_prev[tt - 1, k, 1],
                                     sd = sqrt(Sigma_prev[tt - 1, k, 1, 1]))
          }
        } else if (d > 1) {
          for (k in seq(K)) {
            phi[, k] <- mvtnorm::dmvnorm(y[[tt]],
                                         mean = mu_prev[tt - 1, k, ],
                                         sigma = Sigma_prev[tt - 1, k, , ])
          }
        }
        temp <- t(t(phi) * pi_prev[tt - 1, ])
      }else if (l > 1) {
        if (d == 1) {
          for (k in seq(K)) {
            phi[, k] <- stats::dnorm(y[[tt]],
                                     mean = mu_prev[tt, k, 1],
                                     sd = sqrt(Sigma_prev[tt, k, 1, 1]))
          }
        } else if (d > 1) {
          for (k in seq(K)) {
            phi[, k] <- mvtnorm::dmvnorm(y[[tt]],
                                         mean = mu_prev[tt, k, ],
                                         sigma = Sigma_prev[tt, k, , ])
          }
        }
        temp <- t(t(phi) * pi_prev[tt, ])
      }
      temp_smooth = temp + lap_smooth_const * apply(temp, 1, min)
      
      
      resp[[tt]] <- temp_smooth / rowSums(temp_smooth)
      resp_weighted[[tt]] = diag(biomass[[tt]]) %*% resp[[tt]]
      zest <- resp %>% purrr::map(~ max.col(.x))
      #M-step
      #M-step pi
      resp_sum[[tt]] <- colSums(resp_weighted[[tt]]) %>%
        unlist() %>%
        matrix(ncol = K, byrow = TRUE)
      for (k in seq(K)){
        pi[tt, k] <- (resp_sum [[tt - 1]][, k] + resp_sum [[tt]][, k]) / (rowSums(resp_sum[[tt - 1]]) + rowSums(resp_sum[[tt]]))
        
      }
      #M-step mu
      y_sum[[tt]] <- crossprod(resp_weighted[[tt]], y[[tt]]) %>%
        unlist() %>% 
        array(c(K, d))
      for (k in seq(K)){
        mu[tt, k, ] <- ((resp_sum[[tt - 1]][, k] * mu[tt - 1, k , ]) +  y_sum[[tt]][k, ]) / (resp_sum[[tt - 1]][, k] + resp_sum[[tt]][, k])
        
      }
      
      
      #M-step Sigma
      mat_sum[[tt]] <- array(NA, c(K, d, d))
      yy[[tt]] <- matrix(NA, dim(y[[tt]])[1], d)
      for (k in seq(K)) {
        for(dd in seq(d)) {
          yy [[tt]][, dd] <- (y[[tt]][, dd]- mu[tt, k, dd])
        }
        mat_sum[[tt]][k, , ] <- crossprod(yy[[tt]], yy[[tt]] * resp_weighted[[tt]][, k]) # YY^T * D * YY
      }
      for (k in seq(K)){
        Sigma[tt, k, , ] <- ((resp_sum[[tt - 1]][, k] * Sigma[tt - 1, k, , ]) + mat_sum[[tt]][k, , ]) / (resp_sum[[tt - 1]] [, k] + resp_sum[[tt]] [, k])
      }
      pi_prev <- pi
      mu_prev <- mu
      Sigma_prev <- Sigma 
    }
    
  }
  dimnames(mu) <- list(NULL, paste0("cluster", 1:K), NULL)
  dimnames(Sigma) <- list(NULL, paste0("cluster", 1:K), NULL, NULL)
  dimnames(pi) <- list(NULL, paste0("cluster", 1:K))
  zest <- resp %>% purrr::map(~ max.col(.x))
  fit_init = list(mu = mu, Sigma = Sigma, pi = pi, resp = resp, zest = zest)
  return(fit_init)
}

### flowkernel + CV functions
init_const <- function (y, K, times_to_sample = 50, points_to_sample = 50){
  num_times <- length(y)
  d <- ncol(y[[1]])
  mu <- array(NA, c(num_times, K, d))
  Sigma <- array(NA, c(num_times, K, d, d))
  pi <- matrix(NA, num_times, K)
  # Repeatedly call Mclust until it gives a non-NULL fit:
  init_fit <- NULL
  while (is.null(init_fit)) {
    # subsample data:
    sampled_times <- sample(num_times, times_to_sample, replace=TRUE)
    if (d == 1) {
      sample_data <- y[sampled_times] %>%
        purrr::map(~ .x[sample(nrow(.x), points_to_sample, replace=TRUE)]) %>% 
        unlist()
    }
    else {
      sample_data <- y[sampled_times] %>%
        purrr::map(~ t(.x[sample(nrow(.x), points_to_sample, replace=TRUE), ])) %>% 
        unlist() %>% 
        matrix(ncol = d, byrow = TRUE)
    }
    
    
    if (d == 1) {
      init_fit <- mclust::Mclust(sample_data, G = K, modelNames = "V")
      for (tt in seq(num_times)) {
        mu[tt, , 1] <- init_fit$parameters$mean
        Sigma[tt, , 1, 1] <- init_fit$parameters$variance$sigmasq
        pi[tt, ] <- init_fit$parameters$pro
      }
    } else if (d > 1) {
      init_fit <- mclust::Mclust(sample_data, G = K, modelNames = "VVV")
      if (is.matrix(init_fit$parameters$mean)){
        for (tt in seq(num_times)) {
          mu[tt, ,] <- t(init_fit$parameters$mean)
          pi[tt, ] <- init_fit$parameters$pro
          Sigma[tt, , , ] <- aperm(init_fit$parameters$variance$sigma, c(3,1,2))
        }
      }
    }
  }
  
  #calculate responsibilities
  resp <- calculate_responsibilities(y, mu, Sigma, pi)
  zest <- resp %>% purrr::map(~ max.col(.x))
  list(mu = mu, Sigma = Sigma, pi = pi, resp = resp, zest = zest)
}

#calculate responsibilities, allowing for unevenly spaced times
calculate_responsibilities <- function(y, mu, Sigma, pi){
  resp <- list() # responsibilities gamma[[t]][i, k]
  log_resp <- list() # log of responsibilities
  d <- ncol(y[[1]])
  K <- ncol(mu)
  num_times <- length(y)
  if (d == 1) {
    for (tt in seq(num_times)) {
      log_phi <- matrix(NA, nrow(y[[tt]]), K)
      for (k in seq(K)) {
        log_phi[, k] <- stats::dnorm(y[[tt]],
                                     mean = mu[tt, k, 1],
                                     sd = sqrt(Sigma[tt, k, 1, 1]), log = TRUE)
      }
      log_temp = t(t(log_phi) + log(pi[tt, ]))
      log_resp[[tt]] = log_temp - matrixStats::rowLogSumExps(log_temp)
      resp[[tt]] = exp(log_resp[[tt]])
    }
  } else if (d > 1) {
    for (tt in seq(num_times)) {
      log_phi <- matrix(NA, nrow(y[[tt]]), K)
      for (k in seq(K)) {
        log_phi[, k] <- mvtnorm::dmvnorm(y[[tt]],
                                         mean = mu[tt, k, ],
                                         sigma = Sigma[tt, k, , ], log = TRUE)
      }
      log_temp = t(t(log_phi) + log(pi[tt, ]))
      log_resp[[tt]] = log_temp - matrixStats::rowLogSumExps(log_temp)
      resp[[tt]] = exp(log_resp[[tt]])
    }
  }
  return(resp) 
}

#kernel_em function modified so distances in time are taken from time values, not indices. Also allows for unevenly spaced data points
kernel_em_dates <- function (y, K, hmu, hSigma, hpi, dates = NULL, num_iter = 10, 
                             biomass = default_biomass(y),
                             initial_fit = init_const(y, K, 50, 50)) {
  num_times <- length(y)
  d <- ncol(y[[1]])
  mu <- initial_fit$mu
  Sigma <- initial_fit$Sigma
  pi <- initial_fit$pi
  if (!is.null(dates)){
    numeric_dates <- as.numeric(dates)
    rescaled_dates <- (numeric_dates - min(numeric_dates)) / 3600
  }
  
  for (l in seq(num_iter)) {
    # E-step
    resp <- calculate_responsibilities(y, mu, Sigma, pi)
    resp_weighted <- purrr::map2(biomass, resp, ~ .y * .x)
    
    # M-step
    
    # M-step-pi
    resp_sum <- purrr::map(resp_weighted, ~ colSums(.x)) %>%
      unlist() %>%
      matrix(ncol = K, byrow = TRUE)
    if (is.null(dates)){
      resp_sum_smooth <- apply(
        resp_sum, 2, function(x)
          stats::ksmooth(1:length(x), x, bandwidth = hpi, x.points = 1:length(x))$y
      )
    } else {
      resp_sum_smooth <- apply(
        resp_sum, 2, function(x) 
          stats::ksmooth(rescaled_dates, x, bandwidth = hpi, x.points = rescaled_dates)$y
      )
    }
    
    pi <- resp_sum_smooth / rowSums(resp_sum_smooth)
    
    # M-step-mu
    y_sum <- purrr::map2(resp_weighted, y, ~ crossprod(.x, .y)) %>% 
      unlist() %>% 
      array(c(K, d, num_times)) %>% 
      aperm(c(3,1,2))
    
    
    if (is.null(dates)){
      y_sum_smoothed <- apply(
        y_sum, 2:3, function(x)
          stats::ksmooth(1:length(x), x, bandwidth = hmu, x.points = 1:length(x))$y
      )
      resp_sum_smooth_mu <- apply(
        resp_sum, 2, function(x)
          stats::ksmooth(1:length(x), x, bandwidth = hmu, x.points = 1:length(x))$y
      )
    } else {
      y_sum_smoothed <- apply(
        y_sum, 2:3, function(x) 
          stats::ksmooth(rescaled_dates, x, bandwidth = hmu, x.points = rescaled_dates)$y
      )
      resp_sum_smooth_mu <- apply(
        resp_sum, 2, function(x) 
          stats::ksmooth(rescaled_dates, x, bandwidth = hmu, x.points = rescaled_dates)$y
      )
    }
    
    
    for (j in seq(d)) {
      mu[, , j] <- y_sum_smoothed[, , j] / resp_sum_smooth_mu
    }
    
    # M-step-Sigma
    mat_sum <- array(NA, c(num_times, K, d, d))
    for (tt in seq(num_times)) {
      yy <- matrix(NA, dim(y[[tt]])[1], d)
      for (k in seq(K)) {
        for(dd in seq(d)) {
          yy [,dd] <- (y[[tt]][, dd] - mu[tt, k, dd])
        }
        mat_sum[tt, k, , ] <- crossprod(yy, yy * resp_weighted[[tt]][, k]) 
      }
    }
    if (is.null(dates)){
      mat_sum_smoothed <- apply(
        mat_sum, 2:4, function(x)
          stats::ksmooth(1:length(x), x, bandwidth = hSigma, x.points = 1:length(x))$y
      )
      resp_sum_smooth_Sigma <- apply(
        resp_sum, 2, function(x)
          stats::ksmooth(1:length(x), x, bandwidth = hSigma, x.points = 1:length(x))$y
      )
    } else {
      mat_sum_smoothed <- apply(
        mat_sum, 2:4, function(x)
          stats::ksmooth(rescaled_dates, x, bandwidth = hSigma, x.points = rescaled_dates)$y
      )
      
      resp_sum_smooth_Sigma <- apply(
        resp_sum, 2, function(x) 
          stats::ksmooth(rescaled_dates, x, bandwidth = hSigma, x.points = rescaled_dates)$y
      )
    }
    
    
    for (j in seq(d))
      for (l in seq(d))
        Sigma[, , j, l] <- mat_sum_smoothed[, , j, l] / resp_sum_smooth_Sigma
  }
  
  zest <- resp %>% purrr::map(~ max.col(.x))
  dimnames(mu) <- list(NULL, paste0("cluster", 1:K), NULL)
  dimnames(Sigma) <- list(NULL, paste0("cluster", 1:K), NULL, NULL)
  dimnames(pi) <- list(NULL, paste0("cluster", 1:K))
  
  list(mu = mu, Sigma = Sigma, pi = pi, resp = resp, zest = zest)
}

cv_log_grid_search_parallel <- function(data, 
                                        K, 
                                        hSigma = 10, 
                                        dates = NULL, 
                                        biomass, 
                                        leave_out_every = 5, 
                                        grid_size = 10) {
  start_time <- Sys.time()
  
  # Decide if we are scaling dates or not
  # (If user provides real date/timestamp, scale; if NULL, treat as already scaled)
  if (is.null(dates)) {
    # Just assign 1, 2, 3, ... for each data row
    numeric_dates <- seq_along(data)
    scale_dates <- FALSE
    
    # If not scaling, treat them as hours directly:
    starting_date <- min(numeric_dates)
    middle_date   <- median(numeric_dates)
    end_cv        <- middle_date - starting_date
    
  } else {
    # Convert real calendar times to numeric
    numeric_dates   <- as.numeric(dates)
    scale_dates     <- TRUE
    
    # If scaling, rescale in hours for the "end_cv" logic
    starting_date   <- min(numeric_dates)
    middle_date     <- median(numeric_dates)
    end_cv          <- (middle_date - starting_date) / 3600 
  }
  
  # Generate the sequence for hmu and hpi
  h_values <- exp(seq(log(end_cv), log(1), length = grid_size))
  h_values <- round(h_values)
  h_values <- unique(h_values)
  h_values <- h_values[h_values > 0]
  
  hmu_vals <- h_values
  hpi_vals <- h_values
  
  # Create a data frame of all parameter combinations
  param_grid <- expand.grid(hmu = hmu_vals, hpi = hpi_vals)
  
  # Setup parallel backend
  ncores <- Sys.getenv("SLURM_CPUS_ON_NODE")
  cl <- makeCluster(as.integer(ncores))
  registerDoParallel(cl)
  
  # Export needed variables/functions
  clusterExport(
    cl, 
    varlist = c("kernel_em_cv", "kernel_em_dates", "kernel_em_predict_M",
                "kernel_em_predict", "calculate_responsibilities", "compute_log_like", 
                "data", "K", "hSigma", "biomass", "leave_out_every", "numeric_dates", "scale_dates"),
    envir = environment()
  )
  
  # Load required packages on the workers
  clusterEvalQ(cl, {
    library(flowkernel) 
  })
  
  # Parallelized computation over parameter grid
  cv_scores <- foreach(i = seq_len(nrow(param_grid)), .combine = rbind) %dopar% {
    hmu <- param_grid$hmu[i]
    hpi <- param_grid$hpi[i]
    
    # Pass numeric_dates and scale_dates into kernel_em_cv
    cv_score <- kernel_em_cv(
      y       = data,
      K       = K,
      dates   = numeric_dates,
      hmu     = hmu,
      hSigma  = hSigma,
      hpi     = hpi,
      biomass = biomass,
      leave_out_every = leave_out_every,
      scale_dates     = scale_dates  # <--- new
    )
    
    data.frame(hmu = hmu, hpi = hpi, cv_score = cv_score)
  }
  
  # Stop cluster
  stopCluster(cl)
  
  # Convert cv_scores to matrix form
  results <- acast(cv_scores, hmu ~ hpi, value.var = "cv_score")
  results_df <- as.data.frame(results)
  
  # Re-label columns and rows for clarity
  colnames(results_df) <- paste0("hpi_", colnames(results_df))
  rownames(results_df) <- paste0("hmu_", rownames(results_df))
  
  # Identify best hmu/hpi
  max_pos  <- which(results == max(results, na.rm = TRUE), arr.ind = TRUE)
  best_hmu <- as.numeric(rownames(results)[max_pos[1]])
  best_hpi <- as.numeric(colnames(results)[max_pos[2]])
  
  # Timing
  end_time  <- Sys.time()
  run_time  <- as.numeric(difftime(end_time, start_time, units = "secs"))
  hours     <- floor(run_time / 3600)
  minutes   <- floor((run_time %% 3600) / 60)
  seconds   <- round(run_time %% 60)
  
  cat("Best hmu:", best_hmu, 
      "\nBest hpi:", best_hpi, 
      "\nHighest CV score:", max(results, na.rm = TRUE), 
      "\nTime Elapsed:", run_time, "seconds =>",
      hours, "hours", minutes, "minutes", seconds, "seconds\n")
  
  results_list <- list(
    results      = results_df,
    best_hmu     = best_hmu,
    best_hpi     = best_hpi,
    max_cv_score = max(results, na.rm = TRUE),
    time_elapsed = run_time
  )
  
  return(results_list)
}

#Prediction function

kernel_em_predict <- function(fit, test_dates, train_dates, hmu, hSigma, hpi) {
  num_test <- length(test_dates)
  K <- ncol(fit$pi)
  d <- dim(fit$mu)[3]
  
  numeric_train_dates <- as.numeric(train_dates)
  numeric_test_dates <- as.numeric(test_dates)
  
  # Rescale dates
  min_date <- min(numeric_train_dates)
  rescaled_train_dates <- (numeric_train_dates - min_date) / 3600
  rescaled_test_dates <- (numeric_test_dates - min_date) / 3600
  
  
  # Predict pi 
  pred_pi <- apply(fit$pi, 2, function(x) 
    stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hpi, x.points = rescaled_test_dates)$y
  )
  
  # Predict mu (means) using apply for smoothing over the second and third dimensions
  pred_mu <- apply(fit$mu, 2:3, function(x) 
    stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hmu, x.points = rescaled_test_dates)$y
  )
  
  # Predict Sigma (covariance matrices) using apply over the second, third, and fourth dimensions
  pred_Sigma <- apply(fit$Sigma, 2:4, function(x) 
    stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hSigma, x.points = rescaled_test_dates)$y
  )
  
  list(mu = pred_mu, Sigma = pred_Sigma, pi = pred_pi)
}


kernel_em_predict_M <- function(fit, 
                                test_dates, 
                                train_dates, 
                                train_data, 
                                train_biomass, 
                                hmu, 
                                hSigma, 
                                hpi,
                                scale_dates = TRUE) {
  num_test <- length(test_dates)
  num_train <- length(train_dates)
  K <- ncol(fit$pi)
  d <- dim(fit$mu)[3]
  mu <- array(NA, c(num_test, K, d))
  Sigma <- array(NA, c(num_test, K, d, d))
  pi <- matrix(NA, num_test, K)
  
  # Convert to numeric
  numeric_train_dates <- as.numeric(train_dates)
  numeric_test_dates  <- as.numeric(test_dates)
  
  if (scale_dates) {
    # If we are given real dates, then we do the original rescaling
    min_date <- min(numeric_train_dates)
    rescaled_train_dates <- (numeric_train_dates - min_date) / 3600
    rescaled_test_dates  <- (numeric_test_dates  - min_date) / 3600
  } else {
    # If dates are already in hours or integer indices, treat them directly
    rescaled_train_dates <- numeric_train_dates
    rescaled_test_dates  <- numeric_test_dates
  }
  
  # 1) Predict pi (mixing proportions) ----------------------------------------
  # Weighted responsibilities (weights = biomass)
  resp_weighted <- purrr::map2(train_biomass, fit$resp, ~ .y * .x)
  # Sum over each mixture component across all observations
  resp_sum <- purrr::map(resp_weighted, ~ colSums(.x)) %>%
    unlist() %>%
    matrix(ncol = K, byrow = TRUE)
  
  # Smooth those sums across time
  resp_sum_smooth <- apply(
    resp_sum, 2, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hpi,
                     x.points = rescaled_test_dates)$y
  )
  
  # Normalize to get pi
  pi <- resp_sum_smooth / rowSums(resp_sum_smooth)
  
  # 2) M-step for mu ---------------------------------------------------------
  # Weighted sum of y's for each mixture component
  y_sum <- purrr::map2(resp_weighted, train_data, ~ crossprod(.x, .y)) %>% 
    unlist() %>% 
    array(c(K, d, num_train)) %>% 
    aperm(c(3,1,2))
  
  # Smooth each dimension across time
  y_sum_smoothed <- apply(
    y_sum, 2:3, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hmu,
                     x.points = rescaled_test_dates)$y
  )
  
  # Smooth the sum of responsibilities (again, but with bandwidth hmu)
  resp_sum_smooth_mu <- apply(
    resp_sum, 2, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hmu,
                     x.points = rescaled_test_dates)$y
  )
  
  # Combine smoothed sums to get mu
  for (j in seq(d)) {
    mu[, , j] <- y_sum_smoothed[, , j] / resp_sum_smooth_mu
  }
  
  # 3) M-step for Sigma ------------------------------------------------------
  mat_sum <- array(NA, c(num_train, K, d, d))
  for (tt in seq(num_train)) {
    # Prepare a matrix for each observation's difference from mu
    yy <- matrix(NA, nrow(train_data[[tt]]), d)
    for (k_idx in seq(K)) {
      for (dd in seq(d)) {
        yy[, dd] <- train_data[[tt]][, dd] - fit$mu[tt, k_idx, dd]
      }
      mat_sum[tt, k_idx, , ] <- crossprod(yy, yy * resp_weighted[[tt]][, k_idx])
    }
  }
  
  # Smooth Sigma the same way
  mat_sum_smoothed <- apply(
    mat_sum, 2:4, function(x)
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hSigma,
                     x.points = rescaled_test_dates)$y
  )
  
  resp_sum_smooth_Sigma <- apply(
    resp_sum, 2, function(x) 
      stats::ksmooth(rescaled_train_dates, x, kernel = "normal", bandwidth = hSigma,
                     x.points = rescaled_test_dates)$y
  )
  
  for (j in seq(d)) {
    for (l in seq(d)) {
      Sigma[, , j, l] <- mat_sum_smoothed[, , j, l] / resp_sum_smooth_Sigma
    }
  }
  
  list(mu = mu, Sigma = Sigma, pi = pi)
}

compute_log_like <- function(y, pred) {
  d <- ncol(y[[1]])  # Dimension of the data
  ntimes <- length(y)  # Number of time points
  K <- dim(pred$pi)[2]  # Number of clusters
  log_like <- numeric(ntimes)
  for (tt in seq(ntimes)) {
    likelihood_tt <- 0
    for (k in seq(K)) {
      if (d == 1) {
        # Univariate case
        likelihood_k <- stats::dnorm(y[[tt]], mean = pred$mu[tt, k, 1], 
                                     sd = sqrt(pred$Sigma[tt, k, 1, 1]))
      } else {
        # Multivariate case
        likelihood_k <- mvtnorm::dmvnorm(y[[tt]], mean = pred$mu[tt, k, ], 
                                         sigma = pred$Sigma[tt, k, , ])
      }
      # Weighted likelihood with mixing proportions
      likelihood_tt <- likelihood_tt + pred$pi[tt, k] * likelihood_k
    }
    log_like[tt] <- sum(log(likelihood_tt))
  }
  
  return(mean(log_like))
}


kernel_em_cv <- function(y, 
                         K, 
                         dates, 
                         hmu, 
                         hSigma, 
                         hpi, 
                         biomass, 
                         leave_out_every = 5,
                         scale_dates = TRUE) 
{
  n <- length(y)
  if (leave_out_every <= 0 || leave_out_every >= n) {
    stop("leave_out_every must be a positive integer smaller than the number of data points")
  }
  
  pred <- vector("list", leave_out_every)
  log_like <- vector("list", leave_out_every)
  
  for (i in seq_len(leave_out_every)) {
    test_indices <- seq(i, n, by = leave_out_every)
    train_data   <- y[-test_indices]
    train_date   <- dates[-test_indices]
    train_biomass <- biomass[-test_indices]
    test_data    <- y[test_indices]
    test_date    <- dates[test_indices]
    test_biomass <- biomass[test_indices]
    
    # rescaling vs. not-rescaling
    if (scale_dates){
      fit_train <- kernel_em_dates(y      = train_data,
                                   K      = K,
                                   hmu    = hmu,
                                   hSigma = hSigma,
                                   hpi    = hpi,
                                   dates  = train_date,
                                   biomass = train_biomass)
      
    } else {
      fit_train <- kernel_em_dates(y      = train_data,
                                   K      = K,
                                   hmu    = hmu,
                                   hSigma = hSigma,
                                   hpi    = hpi,
                                   biomass = train_biomass)
    }
    # Pass `scale_dates` here so `kernel_em_predict_M()` knows whether to rescale
    pred[[i]] <- kernel_em_predict_M(fit           = fit_train,
                                     test_dates    = test_date,
                                     train_dates   = train_date,
                                     train_data    = train_data,
                                     train_biomass = train_biomass,
                                     hmu           = hmu,
                                     hSigma        = hSigma,
                                     hpi           = hpi,
                                     scale_dates   = scale_dates)
    
    log_like[[i]] <- compute_log_like(test_data, pred[[i]])
  }
  
  avg_log_likelihood <- mean(unlist(log_like))
  return(avg_log_likelihood)
}
