library(flowkernel)
library(mclust)

## DISAPPEAR SIMULATION ##

d <- 1; K <- 2; ntimes <- 200


disappearance_duration <- 20 # Set the current disappearance duration
transition_duration <- 20

ex1 <- list(
  mu_function = function(times) {
    mu <- array(NA, c(ntimes, K, d))
    mu[, , 1] <- cbind(sin(2 * pi * times / 30), 2)
    mu
  },
  Sigma_function = function(times) {
    Sigma <- array(NA, c(ntimes, K, 1, 1))
    Sigma[, , 1, 1] <- 0.25
    Sigma
  },
  pi_function = function(times) {
    pi2 <- rep(0.5, length(times))
    pi2[times <= transition_duration] <- 0.5 * (1 - times[times <= transition_duration] / transition_duration)
    pi2[(times > transition_duration) & (times <= transition_duration + disappearance_duration)] <- 0
    pi2[times > transition_duration + disappearance_duration] <- 0.5 * 
      (times[times > transition_duration + disappearance_duration] - 
         (transition_duration + disappearance_duration)) / transition_duration
    pi2[pi2 > 0.5] <- 0.5
    
    pi1 <- 1 - pi2
    cbind(pi1, pi2)
  },
  num_points = rep(40, ntimes)
)


# Generate the data for this particular disappearance duration
ex1$dat <- generate_smooth_gauss_mix(ex1$mu_function,
                                     ex1$Sigma_function,
                                     ex1$pi_function,
                                     ex1$num_points)
dis_plot <- plot_data(ex1$dat$y)

# fig is the ggplot object returned by plot_data(...)
dis_plot_clean <- dis_plot +
  labs(title = NULL, x = NULL, y = NULL) +      # drop titles / labels
  theme_minimal() +                             # white background
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15)
  )

ggsave(
  filename = "plot_duration_20.png",
  plot     = dis_plot_clean,
  width    = 1026, height = 673,
  units    = "px", dpi = 100,
  bg       = "white"
)


### Bigger axis labels
## ── enlarge ticks & add “Time” label ──────────────────────────────────────
big_axis_text  <- 18   # tick labels
big_axis_title <- 18   # “Time” axis title

dis_plot_clean <- dis_plot +
  labs(x = "Time", y = NULL, title = NULL) +      # x-label added
  theme_minimal() +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
    axis.text        = element_text(size = big_axis_text),
    axis.title       = element_text(size = big_axis_title, face = "plain")
  )

ggsave(
  filename = "dis_plot_duration_60.png",
  plot     = dis_plot_clean,
  width    = 1026, height = 673,
  units    = "px", dpi = 100,
  bg       = "white"
)

## End bigger axis labels


constant_fit  <- const_fit(ex1$dat$y, K = 2)
const_plot    <- plot_data_and_model(ex1$dat$y, constant_fit$zest,
                                     mu = constant_fit$mu)

clean_and_save(const_plot, "dis_20_const_fit.png",
               mean_width = 1)   # e.g. make mean lines noticeably thicker


## Kernel-em fit
flowkernel::cv_grid_search(y = ex1$dat$y, K = 2) # search for params
set.seed(1243)
kernel_em_fit <- kernel_em(ex1$dat$y, K = 2, hmu = 11, hSigma = 15, hpi = 12, num_iter = 10)
plot_data_and_model(ex1$dat$y, kernel_em_fit$zest, mu = kernel_em_fit$mu)
kernel_em_fit_plot <- plot_data_and_model(ex1$dat$y, kernel_em_fit$zest, mu = kernel_em_fit$mu)

clean_and_save(kernel_em_fit_plot, "dis_20_kernel_fit.png",
               mean_width = 1)  # e.g. make mean lines noticeably thicker

## Hungarian fit

hungarian_fit_dis <- hungarian_fit(y = ex1$dat$y, K = 2)
hung_plot <-plot_data_and_model(ex1$dat$y, hungarian_fit_dis$zest, mu = hungarian_fit_dis$mu)

clean_and_save(hung_plot, "dis_20_hung_fit.png",
               mean_width = 1)   # e.g. make mean lines noticeably thicker

## Clean and save function
clean_and_save <- function(p, file) {
  p_clean <- p +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15)
    )
  
  ggsave(
    filename = file,
    plot     = p_clean,
    width    = 1026, height = 673,
    units    = "px", dpi = 100, bg = "white"
  )
}

### Clean and save bigger axis titles

## ── overwrite helper ──────────────────────────────────────────────────────
# ── helper: make ALL line / hline layers thicker ──────────────────────────
thicken_means <- function(p, new_width = 2) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, c("GeomLine", "GeomPath", "GeomHline"))) {
      
      # ensure the lists exist
      if (is.null(p$layers[[i]]$aes_params))  p$layers[[i]]$aes_params  <- list()
      if (is.null(p$layers[[i]]$geom_params)) p$layers[[i]]$geom_params <- list()
      
      # modern slot (ggplot2 ≥3.4)
      p$layers[[i]]$aes_params$linewidth  <- new_width
      p$layers[[i]]$geom_params$linewidth <- new_width
      
      # legacy / fallback slot (ggplot2 <3.4)
      p$layers[[i]]$aes_params$size  <- new_width
      p$layers[[i]]$geom_params$size <- new_width
    }
  }
  p
}

# ── clean-and-save with big axis text AND thicker means ───────────────────
clean_and_save <- function(p, file,
                           axis_text_size  = 18,
                           axis_title_size = 18,
                           mean_width      = 3) {
  
  p_clean <- p |>
    thicken_means(mean_width) +                 # << fatten mean lines
    labs(x = "Time", y = NULL, title = NULL) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
      axis.text        = element_text(size = axis_text_size),
      axis.title       = element_text(size = axis_title_size, face = "plain")
    )
  
  ggsave(
    filename = file,
    plot     = p_clean,
    width    = 1026, height = 673,
    units    = "px", dpi = 100, bg = "white"
  )
  
  invisible(p_clean)
}



### Fitting Functions ###

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

### Cluster intersection
library(flowkernel)
library(ggplot2)
library(purrr)
library(tibble)
library(dplyr)

## ------------------------------------------------------------------
## 1.  Choose the intersection level you want to visualise
## ------------------------------------------------------------------
intersect_cluster <- 2.9            # ← change this number as needed

## ------------------------------------------------------------------
## 2.  Build the functions for this model
## ------------------------------------------------------------------
ntimes <- 200;  d <- 1;  K <- 2

mu_function <- function(times) {
  mu <- array(NA_real_, c(ntimes, K, d))
  mu[, , 1] <- cbind(
    cos(2 * pi * times / 30) + intersect_cluster,          # cluster 1
    atan(2 * pi * times / 30) + 2.2                        # cluster 2
  )
  mu
}

Sigma_function <- function(times) {
  array(0.20, c(ntimes, K, d, d))                          # constant σ² = 0.20
}

pi_function <- function(times) {
  pi1 <- seq(0.2, 0.8, length.out = length(times))
  cbind(pi1, 1 - pi1)
}

num_points <- rep(40, ntimes)

## ------------------------------------------------------------------
## 3.  Generate the data and grab the raw scatter plot
## ------------------------------------------------------------------
dat <- generate_smooth_gauss_mix(
  mu_function,
  Sigma_function,
  pi_function,
  num_points
)

raw_plot <- plot_data(dat$y)

## ------------------------------------------------------------------
## 4.  Clean styling: no titles / labels, white bg with light grid
## ------------------------------------------------------------------
plot_clean <- raw_plot +
  labs(title = NULL, x = "Time", y = NULL) +     # x-axis label
  theme_minimal() +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
    panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
    axis.text        = element_text(size = 18),  # tick numbers
    axis.title       = element_text(size = 18, face = "plain")  # “Time”
  )


## ------------------------------------------------------------------
## 5.  Save at the same dimensions you’ve been using
## ------------------------------------------------------------------
ggsave(
  filename = sprintf("plot_intersection_%s.png", intersect_cluster),
  plot     = plot_clean,
  width    = 1026, height = 673,
  units    = "px", dpi = 100,
  bg       = "white"
)

## Loop ##
library(flowkernel)
library(ggplot2)
library(purrr)
library(tibble)
library(dplyr)

## -------------------------------------------------------------
##  global settings (same for every intersection level)
## -------------------------------------------------------------
ntimes      <- 200
d           <- 1
K           <- 2
times_index <- 1:ntimes                      # the “times” each function sees
num_points  <- rep(40, ntimes)

Sigma_function <- function(t) {
  array(0.20, c(length(t), K, d, d))         # constant variance
}

pi_function <- function(t) {
  pi1 <- seq(0.2, 0.8, length.out = length(t))
  cbind(pi1, 1 - pi1)
}

## -------------------------------------------------------------
##  helper: clean-up & save a ggplot
## -------------------------------------------------------------
save_clean_plot <- function(p, filename) {
  p_clean <- p +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      legend.position  = "none",
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15)
    )
  
  ggsave(
    filename = filename,
    plot     = p_clean,
    width    = 1026, height = 673,
    units    = "px", dpi = 100,
    bg       = "white"
  )
  message("saved ", filename)
}

## -------------------------------------------------------------
##  loop over intersection levels 0.5, 0.8, …, 2.9
## -------------------------------------------------------------
for (icl in seq(0.5, 2.9, by = 0.3)) {
  
  ## intersection-specific mean function ------------------------
  mu_function <- function(t) {
    mu <- array(NA_real_, c(length(t), K, d))
    mu[, , 1] <- cbind(
      cos(2 * pi * t / 30) + icl,            # cluster 1 mean
      atan(2 * pi * t / 30) + 2.2            # cluster 2 mean
    )
    mu
  }
  
  ## generate one dataset --------------------------------------
  dat <- generate_smooth_gauss_mix(
    mu_function,
    Sigma_function,
    pi_function,
    num_points
  )
  
  ## make raw scatter plot -------------------------------------
  p_raw <- plot_data(dat$y)
  
  ## clean + save ----------------------------------------------
  fn <- sprintf("plot_intersection_%s.png",
                gsub("\\.", "p", formatC(icl, format = "f", digits = 1)))
  save_clean_plot(p_raw, fn)
}



