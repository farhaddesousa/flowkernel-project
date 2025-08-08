## @knitr models

# model_functions.R
library(simulator)

# Define the model-generating function for a disappearing cluster
make_disappearance_model <- function(transition_duration, disappearance_duration) {
  ntimes <- 200
  d <- 1
  K <- 2
  
  mu_function <- function(times) {
    mu <- array(NA, c(ntimes, K, d))
    mu[, , 1] <- cbind(sin(2 * pi * times / 30), 2)
    mu
  }
  
  Sigma_function <- function(times) {
    Sigma <- array(NA, c(ntimes, K, 1, 1))
    Sigma[, , 1, 1] <- 0.25
    Sigma
  }
  
  pi_function <- function(times) {
    pi2 <- rep(0.5, length(times))
    pi2[times <= transition_duration] <- 0.5 * (1 - times[times <= transition_duration] / transition_duration)
    pi2[(times > transition_duration) & (times <= transition_duration + disappearance_duration)] <- 0
    pi2[times > transition_duration + disappearance_duration] <- 0.5 * 
      (times[times > transition_duration + disappearance_duration] - 
         (transition_duration + disappearance_duration)) / transition_duration
    pi2[pi2 > 0.5] <- 0.5
    
    pi1 <- 1 - pi2
    cbind(pi1, pi2)
  }
  
  num_points <- rep(40, ntimes)
  
  # Create the model object
  new_model(name = sprintf("disappearance_model_dur_%s", disappearance_duration),
            label = sprintf("Disappearance Duration = %s", disappearance_duration),
            params = list(mu_function = mu_function,
                          Sigma_function = Sigma_function,
                          pi_function = pi_function,
                          num_points = num_points,
                          transition_duration = transition_duration,
                          disappearance_duration = disappearance_duration),
            simulate = function(mu_function, Sigma_function, pi_function, num_points, nsim) {
              # Generate data
              replicate(nsim, {
                dat <- generate_smooth_gauss_mix(mu_function,
                                                 Sigma_function,
                                                 pi_function,
                                                 num_points)
                list(y = dat$y, dates = dat$dates, biomass = dat$biomass, z = dat$z)
              }, simplify = FALSE)
            })
}

#Make 3-d disappearance model:

make_disappearance_model_3d <- function(transition_duration, disappearance_duration, disappearing_cluster = 2) {
  ntimes <- 200
  d <- 3
  K <- 4
  
  mu_function <- function(times) {
    mu <- array(NA, c(ntimes, K, d))
    mu[, , 1] <- cbind(
      0.5 * cos(2 * pi * times / 30),
      0.3 * sin(2 * pi * times / 30),
      sin(2 * pi * times / 30),
      -3
    )
    mu[, , 2] <- cbind(
      0.3 * sin(2 * pi * times / 30),
      2,
      -1,
      0.6 * cos(2 * pi * times / 30)
    )
    mu[, , 3] <- cbind(
      2,
      0.7 * cos(2 * pi * times / 30),
      0.4 * sin(2 * pi * times / 30),
      1
    )
    mu
  }
  
  Sigma_function <- function(times) {
    Sigma <- array(0, c(ntimes, K, d, d))
    for (k in 1:K) {
      for (t in 1:ntimes) {
        Sigma[t, k, , ] <- diag(rep(0.1, d))
      }
    }
    Sigma
  }
  
  pi_function <- function(times) {
    pi_matrix <- matrix(NA, nrow = length(times), ncol = K)
    pi1 <- seq(0.2, 0.3, length = length(times))
    pi_matrix[, 1] <- pi1
    pi_matrix[, 2] <- pi1
    pi_matrix[, 3] <- (2 * pi1) / 3
    pi_matrix[, 4] <- 1 - rowSums(pi_matrix[, 1:3])
    
    # Adjust the disappearing cluster's probabilities
    pi_disappear <- pi_matrix[, disappearing_cluster]
    pi_disappear_initial <- pi_disappear
    
    t_decrease_start <- transition_duration
    t_decrease_end <- transition_duration + disappearance_duration / 2
    t_increase_start <- t_decrease_end
    t_increase_end <- transition_duration + disappearance_duration
    
    # Decrease to zero
    indices_decrease <- times >= t_decrease_start & times <= t_decrease_end
    pi_disappear[indices_decrease] <- pi_disappear_initial[indices_decrease] * (
      1 - (times[indices_decrease] - t_decrease_start) / (t_decrease_end - t_decrease_start)
    )
    
    # Stay zero
    indices_zero <- times > t_decrease_end & times < t_increase_start
    pi_disappear[indices_zero] <- 0
    
    # Increase back to initial value
    indices_increase <- times >= t_increase_start & times <= t_increase_end
    pi_disappear[indices_increase] <- pi_disappear_initial[indices_increase] * (
      (times[indices_increase] - t_increase_start) / (t_increase_end - t_increase_start)
    )
    
    # Update the disappearing cluster's probabilities
    pi_matrix[, disappearing_cluster] <- pi_disappear
    
    # Recalculate the last cluster's probabilities
    pi_matrix[, 4] <- 1 - rowSums(pi_matrix[, 1:3])
    
    # Ensure no negative probabilities
    pi_matrix[pi_matrix < 0] <- 0
    
    # Normalize the probabilities
    row_sums <- rowSums(pi_matrix)
    pi_matrix <- pi_matrix / row_sums
    
    return(pi_matrix)
  }
  
  num_points <- rep(150, ntimes)
  
  # Create the model object
  new_model(
    name = sprintf("disappearance_model_3d_cluster_%s", disappearing_cluster),
    label = sprintf("Disappearance of Cluster %s", disappearing_cluster),
    params = list(
      mu_function = mu_function,
      Sigma_function = Sigma_function,
      pi_function = pi_function,
      num_points = num_points,
      transition_duration = transition_duration,
      disappearance_duration = disappearance_duration,
      disappearing_cluster = disappearing_cluster
    ),
    simulate = function(mu_function, Sigma_function, pi_function, num_points, nsim) {
      # Generate data
      replicate(nsim, {
        dat <- generate_smooth_gauss_mix(
          mu_function,
          Sigma_function,
          pi_function,
          num_points
        )
        list(y = dat$y, dates = dat$dates, biomass = dat$biomass, z = dat$z)
      }, simplify = FALSE)
    }
  )
}




# Define the model-generating function for intersecting clusters
make_intersection_model <- function(intersect_cluster) {
  ntimes <- 200
  d <- 1
  K <- 2
  
  # Define the functions to generate parameters
  mu_function <- function(times) {
    mu <- array(NA, c(ntimes, K, d))
    # Define means to vary with intersect_cluster value
    mu[, , 1] <- cbind(cos(2 * pi * times / 30) + intersect_cluster, atan(2 * pi * times / 30) + 2.2)
    mu
  }
  
  Sigma_function <- function(times) {
    Sigma <- array(NA, c(ntimes, K, 1, 1))
    Sigma[, , 1, 1] <- 0.20
    Sigma
  }
  
  pi_function <- function(times) {
    pi1 <- seq(0.2, 0.8, length = length(times))
    pi2 <- 1 - pi1
    cbind(pi1, pi2)
  }
  
  num_points <- rep(40, ntimes)
  
  # Create the model object
  new_model(name = sprintf("intersection_model_intersect_%s", intersect_cluster),
            label = sprintf("Intersection Cluster = %s", intersect_cluster),
            params = list(mu_function = mu_function,
                          Sigma_function = Sigma_function,
                          pi_function = pi_function,
                          num_points = num_points,
                          intersect_cluster = intersect_cluster),
            simulate = function(mu_function, Sigma_function, pi_function, num_points, nsim) {
              # Generate data
              replicate(nsim, {
                dat <- generate_smooth_gauss_mix(mu_function,
                                                 Sigma_function,
                                                 pi_function,
                                                 num_points)
                list(y = dat$y, dates = dat$dates, biomass = dat$biomass, z = dat$z)
              }, simplify = FALSE)
            })
}

#New intersect model - flip
make_intersection_model_flip <- function(intersect_cluster) {
  ntimes <- 200
  d <- 1
  K <- 2
  period <- 50
  
  # Define cluster means that cross as intersect_cluster increases.
  # Cluster 1: a sine wave.
  # Cluster 2: a negative sine wave shifted by an offset that decreases with intersect_cluster.
  mu_function <- function(times) {
    mu <- array(NA, c(ntimes, K, d))
    # Cluster 1: sine wave
    mu[, 1, 1] <- sin(2 * pi * times / period)
    # Cluster 2: negative sine wave plus an offset (the offset decreases as intersect_cluster increases)
    offset <- 3 - intersect_cluster  # e.g., if intersect_cluster = 0.5 then offset = 2.5; if =3 then offset = 0
    mu[, 2, 1] <- -sin(2 * pi * times / period) + offset
    mu
  }
  
  # Use a slightly larger variance to further increase overlap when means come close.
  Sigma_function <- function(times) {
    Sigma <- array(NA, c(ntimes, K, 1, 1))
    Sigma[, , 1, 1] <- 0.3
    Sigma
  }
  
  # Set the mixing proportions constant: first cluster = 0.7, second cluster = 0.3
  pi_function <- function(times) {
    cbind(rep(0.7, length(times)), rep(0.3, length(times)))
  }
  
  num_points <- rep(20, ntimes)
  
  new_model(name = sprintf("intersection_model_intersect_%s", intersect_cluster),
            label = sprintf("Intersection Cluster = %s", intersect_cluster),
            params = list(mu_function = mu_function,
                          Sigma_function = Sigma_function,
                          pi_function = pi_function,
                          num_points = num_points,
                          intersect_cluster = intersect_cluster),
            simulate = function(mu_function, Sigma_function, pi_function, num_points, nsim) {
              replicate(nsim, {
                dat <- generate_smooth_gauss_mix(mu_function,
                                                 Sigma_function,
                                                 pi_function,
                                                 num_points)
                list(y = dat$y, dates = dat$dates, biomass = dat$biomass, z = dat$z)
              }, simplify = FALSE)
            })
}





generate_smooth_gauss_mix <- function(mu_function,
                                      Sigma_function,
                                      pi_function,
                                      num_points,
                                      start_date = "2017-05-28 21:00:00 UTC") {
  times <- seq_along(num_points)
  mu <- mu_function(times)
  Sigma <- Sigma_function(times)
  pi <- pi_function(times)
  K <- ncol(pi) # number of components
  d <- dim(mu)[3]
  dimnames(mu) <- list(NULL, paste0("cluster", 1:K), NULL)
  
  z <- list() # z[[t]][i] = class of point i at time t
  y <- list() # y[[t]][i,] = d-vector of point i at time t
  biomass <- list() # biomass[[t]] = biomass of particles in each bin at time t
  for (t in times) {
    z[[t]] <- apply(stats::rmultinom(num_points[t], 1, pi[t, ]) == 1, 2, which)
    y[[t]] <- matrix(NA, num_points[t], d)
    biomass[[t]] <- abs(rnorm(num_points[t], mean = 0.01, sd = sqrt(0.001)))
    for (k in 1:K) {
      ii <- z[[t]] == k # index of points in component k at time t
      if (sum(ii) == 0) next
      if (d == 1)
        y[[t]][ii, ] <- stats::rnorm(n = sum(ii),
                                     mean = mu[t, k, ],
                                     sd = Sigma[t, k, , ])
      else
        y[[t]][ii, ] <- mvtnorm::rmvnorm(n = sum(ii),
                                         mean = mu[t, k, ],
                                         sigma = Sigma[t, k, , ])
    }
  }
  
  # Generate dates spaced 1 hour apart
  start_time <- as.POSIXct(start_date, tz = "UTC")
  dates <- seq(start_time, by = "hour", length.out = length(y))
  
  list(y = y, z = z, mu = mu, Sigma = Sigma, pi = pi, biomass = biomass, dates = dates)
}
