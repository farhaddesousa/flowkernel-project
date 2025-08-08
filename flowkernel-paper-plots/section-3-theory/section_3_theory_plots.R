library(mvtnorm)
library(ggplot2)
library(tibble)
library(dplyr)
library(purrr)

# 1) Data-generation function:
#    Note: Sigma_function must return a 4D array of dimensions
#    [length(times), K, d, d] even for univariate data.
generate_smooth_gauss_mix_cross <- function(mu_function,
                                            Sigma_function,
                                            pi_function,
                                            times,
                                            num_points) {
  # Evaluate functions at all time points:
  mu    <- mu_function(times)     # dimensions: [length(times), K, d]
  Sigma <- Sigma_function(times)    # dimensions: [length(times), K, d, d]
  pi    <- pi_function(times)       # dimensions: [length(times), K]
  
  K <- ncol(pi)        # Number of mixture components
  d <- dim(mu)[3]      # Dimensionality of the data
  N <- length(times)   # Number of time points
  
  # Pre-allocate lists:
  z       <- vector("list", N)
  y       <- vector("list", N)
  biomass <- vector("list", N)
  
  # Loop over time points:
  for (i in seq_len(N)) {
    n_i <- num_points[i]
    
    # Draw component labels based on mixture probabilities:
    z[[i]] <- apply(stats::rmultinom(n_i, 1, pi[i, ]) == 1, 2, which)
    
    # Allocate storage for the data at this time:
    y[[i]] <- matrix(NA, n_i, d)
    biomass[[i]] <- abs(rnorm(n_i, mean = 0.01, sd = sqrt(0.001)))
    
    # Generate data for each mixture component:
    for (k in seq_len(K)) {
      these <- z[[i]] == k
      if (!any(these)) next
      
      if (d == 1) {
        # Univariate normal: interpret Sigma[i,k,1,1] as the variance.
        y[[i]][these, ] <- stats::rnorm(
          n    = sum(these),
          mean = mu[i, k, ],
          sd   = sqrt(Sigma[i, k, 1, 1])
        )
      } else {
        # Multivariate normal:
        y[[i]][these, ] <- mvtnorm::rmvnorm(
          n     = sum(these),
          mean  = mu[i, k, ],
          sigma = Sigma[i, k, , ]
        )
      }
    }
  }
  
  list(
    y       = y,
    z       = z,
    mu      = mu,
    Sigma   = Sigma,
    pi      = pi,
    biomass = biomass,
    times   = times
  )
}

# 2) Plotting function:
#    For the univariate case, we construct a tibble with an explicit “V1” column.
plot_data_cross <- function(dat) {
  if (!("times" %in% names(dat))) {
    stop("dat must include a 'times' element.")
  }
  
  the_times <- dat$times
  the_y     <- dat$y
  
  # For univariate data, each y[[i]] is a matrix with one column.
  df <- map2_dfr(
    .x = the_y,
    .y = the_times,
    .f = function(mat, tval) {
      tibble(V1 = as.vector(mat[, 1]), time = tval)
    }
  )
  
  ggplot(df, aes(x = time, y = V1)) +
    geom_point(alpha = 0.3) +
    labs(
      x = "Time",
      y = "V1"
      #title = "Raw Data (Mixture of Gaussians Over Time)"
    ) +
    scale_x_continuous(limits = c(-1, 1))
}

# 3) Example usage:
# Define times from -1 to 1 in steps of 0.01.
times <- seq(-1, 1, by = 0.01)
ntimes <- length(times)
num_points <- rep(40, ntimes)  # 40 points per time step

# Set number of clusters and dimension.
K <- 2  # two clusters
d <- 1  # univariate data

# Define the functions.
mu_function <- function(t) {
  # Create an array of dimensions [length(t), K, d].
  arr <- array(NA_real_, dim = c(length(t), K, d))
  arr[, 1, 1] <-  t   # For cluster 1, mean equals time.
  arr[, 2, 1] <- -t   # For cluster 2, mean equals negative time.
  arr
}

Sigma_function <- function(t) {
  # Return a constant variance of 0.1 for each component.
  array(0.01, dim = c(length(t), K, d, d))
}

pi_function <- function(t) {
  # Both clusters have constant proportions (0.5 each).
  matrix(0.5, nrow = length(t), ncol = K)
}

# Generate the data.
ex_data <- generate_smooth_gauss_mix_cross(
  mu_function    = mu_function,
  Sigma_function = Sigma_function,
  pi_function    = pi_function,
  times          = times,
  num_points     = num_points
)

plot_data_cross <- function(dat) {
  if (!("times" %in% names(dat))) {
    stop("dat must include a 'times' element.")
  }
  
  the_times <- dat$times
  the_y     <- dat$y
  
  # For univariate data, each y[[i]] is a matrix with one column.
  df <- map2_dfr(
    .x = the_y,
    .y = the_times,
    .f = function(mat, tval) {
      tibble(V1 = as.vector(mat[, 1]), time = tval)
    }
  )
  
  ggplot(df, aes(x = time, y = V1)) +
    geom_point(color = "black", alpha = 0.3) +  # Match aesthetics
    labs(
      x = expression("Time"),
      y = expression("V1")
    ) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +  # Match the y-axis range
    theme_minimal()
}

# Plot the generated data.
plot_data_cross(ex_data)


library(ggplot2)

p_raw <- plot_data_cross(ex_data) + labs(y = NULL)

ggsave(
  "raw_crossing_ex.png",   # overwrite or rename as needed
  plot   = p_raw,
  width  = 1026, height = 673,  # identical pixel size to the old file
  units  = "px",
  dpi    = 100,                 # bump to 600 for print-ready output
  bg     = "white"
)



### Bigger numbers on axis ###
library(ggplot2)
library(purrr)
library(tibble)

plot_data_cross <- function(dat,
                            axis_text_size  = 18,   # << bigger tick numbers
                            axis_title_size = 18) { # << bigger “Time” label
  if (!("times" %in% names(dat))) {
    stop("dat must include a 'times' element.")
  }
  
  df <- map2_dfr(dat$y, dat$times,
                 ~ tibble(V1 = as.vector(.x[, 1]), time = .y))
  
  ggplot(df, aes(time, V1)) +
    geom_point(colour = "black", alpha = 0.3) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(x = "Time", y = NULL) +             # << x-axis label
    theme_minimal() +
    theme(
      axis.text  = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size, face = "plain")
    )
}

p_raw <- plot_data_cross(ex_data)

ggsave(
  "raw_crossing_theory.png",
  plot   = p_raw,
  width  = 1026, height = 673,  # keep identical pixel size
  units  = "px",
  dpi    = 100,                 # unchanged
  bg     = "white"
)


### End bigger numbers ### 

######################################################################
# 1) Data-generation function:
#    Note: Sigma_function must return a 4D array of dimensions
#    [length(times), K, d, d] even for univariate data.
generate_smooth_gauss_mix_cross <- function(mu_function,
                                            Sigma_function,
                                            pi_function,
                                            times,
                                            num_points) {
  # Evaluate functions at all time points:
  mu    <- mu_function(times)     # dimensions: [length(times), K, d]
  Sigma <- Sigma_function(times)    # dimensions: [length(times), K, d, d]
  pi    <- pi_function(times)       # dimensions: [length(times), K]
  
  K <- ncol(pi)        # Number of mixture components
  d <- dim(mu)[3]      # Dimensionality of the data
  N <- length(times)   # Number of time points
  
  # Pre-allocate lists:
  z       <- vector("list", N)
  y       <- vector("list", N)
  biomass <- vector("list", N)
  
  # Loop over time points:
  for (i in seq_len(N)) {
    n_i <- num_points[i]
    
    # Draw component labels based on mixture probabilities:
    z[[i]] <- apply(stats::rmultinom(n_i, 1, pi[i, ]) == 1, 2, which)
    
    # Allocate storage for the data at this time:
    y[[i]] <- matrix(NA, n_i, d)
    biomass[[i]] <- abs(rnorm(n_i, mean = 0.01, sd = sqrt(0.001)))
    
    # Generate data for each mixture component:
    for (k in seq_len(K)) {
      these <- z[[i]] == k
      if (!any(these)) next
      
      if (d == 1) {
        # Univariate normal: interpret Sigma[i,k,1,1] as the variance.
        y[[i]][these, ] <- stats::rnorm(
          n    = sum(these),
          mean = mu[i, k, ],
          sd   = sqrt(Sigma[i, k, 1, 1])
        )
      } else {
        # Multivariate normal:
        y[[i]][these, ] <- mvtnorm::rmvnorm(
          n     = sum(these),
          mean  = mu[i, k, ],
          sigma = Sigma[i, k, , ]
        )
      }
    }
  }
  
  list(
    y       = y,
    z       = z,
    mu      = mu,
    Sigma   = Sigma,
    pi      = pi,
    biomass = biomass,
    times   = times
  )
}

# 2) Plotting function:
#    For the univariate case, we construct a tibble with an explicit “V1” column.
plot_data_cross <- function(dat) {
  if (!("times" %in% names(dat))) {
    stop("dat must include a 'times' element.")
  }
  
  the_times <- dat$times
  the_y     <- dat$y
  
  # For univariate data, each y[[i]] is a matrix with one column.
  df <- map2_dfr(
    .x = the_y,
    .y = the_times,
    .f = function(mat, tval) {
      tibble(V1 = as.vector(mat[, 1]), time = tval)
    }
  )
  
  ggplot(df, aes(x = time, y = V1)) +
    geom_point(alpha = 0.3) +
    labs(
      x = "Time",
      y = "V1"
      #title = "Raw Data (Mixture of Gaussians Over Time)"
    ) +
    scale_x_continuous(limits = c(-1, 1))
}

# 3) Example usage:
# Define times from -1 to 1 in steps of 0.01.
times <- seq(-1, 1, by = 0.01)
ntimes <- length(times)
num_points <- rep(40, ntimes)  # 40 points per time step

# Set number of clusters and dimension.
K <- 2  # two clusters
d <- 1  # univariate data

# Define the functions.
mu_function <- function(t) {
  # Create an array of dimensions [length(t), K, d].
  arr <- array(NA_real_, dim = c(length(t), K, d))
  arr[, 1, 1] <-  t   # For cluster 1, mean equals time.
  arr[, 2, 1] <- -t   # For cluster 2, mean equals negative time.
  arr
}

Sigma_function <- function(t) {
  # Return a constant variance of 0.1 for each component.
  array(0.01, dim = c(length(t), K, d, d))
}

pi_function <- function(t) {
  # Both clusters have constant proportions (0.5 each).
  matrix(0.5, nrow = length(t), ncol = K)
}

# Generate the data.
ex_data <- generate_smooth_gauss_mix_cross(
  mu_function    = mu_function,
  Sigma_function = Sigma_function,
  pi_function    = pi_function,
  times          = times,
  num_points     = num_points
)
plot_data_cross <- function(dat) {
  if (!("times" %in% names(dat))) {
    stop("dat must include a 'times' element.")
  }
  
  the_times <- dat$times
  the_y     <- dat$y
  
  # For univariate data, each y[[i]] is a matrix with one column.
  df <- map2_dfr(
    .x = the_y,
    .y = the_times,
    .f = function(mat, tval) {
      tibble(V1 = as.vector(mat[, 1]), time = tval)
    }
  )
  
  ggplot(df, aes(x = time, y = V1)) +
    geom_point(color = "black", alpha = 0.3) +  # Match aesthetics
    labs(
      x = expression("Time"),
      y = expression("V1")
    ) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +  # Match the y-axis range
    theme_minimal()
}

# Plot the generated data.
plot_data_cross(ex_data)


### 
library(ggplot2)

# Create a sequence of x values from -1 to 1
x_vals <- seq(-1, 1, length.out = 100)

# Define the y values for both plots
df1 <- data.frame(x = x_vals, y1 = x_vals, y2 = -x_vals)  # Linear functions
df2 <- data.frame(x = x_vals, y1 = abs(x_vals), y2 = -abs(x_vals))  # Absolute value functions

# First plot: Linear functions crossing
p1 <- ggplot(df1, aes(x = x)) +
  geom_line(aes(y = y1), color = "red", size = 1.2) +       # Red solid line
  geom_line(aes(y = y2), color = "blue", linetype = "dashed", size = 1.2) + # Blue dashed line
  xlim(-1, 1) + ylim(-1, 1) +
  labs(
    x = expression("Time"), 
    y = expression("Cluster Means")
  ) +
  theme_minimal()


p1
# Save first plot
ggsave("linear_functions.png", p1, width = 6, height = 4, dpi = 300)

# Second plot: Absolute value functions
p2 <- ggplot(df2, aes(x = x)) +
  geom_line(aes(y = y1), color = "blue", linetype = "dashed", size = 1.2) +  # Blue dashed abs(x)
  geom_line(aes(y = y2), color = "red", size = 1.2) +      # Red solid -abs(x)
  xlim(-1, 1) + ylim(-1, 1) +
  labs(
    x = expression("Time"), 
    y = expression("Cluster Means")
  ) +
  theme_minimal()
p2

library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

#### 1) Define the time grid and mixture parameters ####
times <- seq(-1, 1, by = 0.01)
ntimes <- length(times)
num_points <- rep(40, ntimes)  # 40 points per time step

K <- 2  # two clusters
d <- 1  # univariate data

Sigma_function <- function(t) {
  # Constant variance of 0.01 for each cluster and time point
  array(0.01, dim = c(length(t), K, d, d))
}

pi_function <- function(t) {
  # Equal mixture weights (0.5, 0.5)
  matrix(0.5, nrow = length(t), ncol = K)
}

#### 2) Define two different mean functions ####
#   a) Linear crossing: mu1 = t, mu2 = -t
mu_function_linear <- function(t) {
  arr <- array(NA_real_, dim = c(length(t), K, d))
  arr[, 1, 1] <- t      # cluster 1: y = t
  arr[, 2, 1] <- -t     # cluster 2: y = -t
  arr
}

#   b) Absolute value crossing: mu1 = |t|, mu2 = -|t|
mu_function_abs <- function(t) {
  arr <- array(NA_real_, dim = c(length(t), K, d))
  arr[, 1, 1] <- abs(t)   # cluster 1: y = |t|
  arr[, 2, 1] <- -abs(t)  # cluster 2: y = -|t|
  arr
}

#### 3) Generate data from each mean specification ####
ex_data_linear <- generate_smooth_gauss_mix_cross(
  mu_function    = mu_function_linear,
  Sigma_function = Sigma_function,
  pi_function    = pi_function,
  times          = times,
  num_points     = num_points
)

ex_data_abs <- generate_smooth_gauss_mix_cross(
  mu_function    = mu_function_abs,
  Sigma_function = Sigma_function,
  pi_function    = pi_function,
  times          = times,
  num_points     = num_points
)

#### 4) Convert each data set into a single tibble with cluster labels ####
# Note that ex_data$z[[i]] is the cluster label (1 or 2) for each point at time[i].
to_df <- function(dat) {
  map2_dfr(
    .x = dat$y,
    .y = seq_along(dat$y),
    .f = function(yi, i) {
      tibble(
        time    = dat$times[i],
        V1      = yi[, 1],
        cluster = factor(dat$z[[i]])  # cluster 1 or 2
      )
    }
  )
}

df_linear <- to_df(ex_data_linear)
df_abs    <- to_df(ex_data_abs)
p_linear <- ggplot(df_linear, aes(x = time, y = V1, color = cluster)) +
  geom_point(alpha = 0.4) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Time",
    y = "V1",
    title = "Linear Means (y = t, y = -t)"
  ) +
  theme_minimal()

p_abs <- ggplot(df_abs, aes(x = time, y = V1, color = cluster)) +
  geom_point(alpha = 0.4) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    x = "Time",
    y = "V1",
    title = "Absolute Value Means (y = |t|, y = -|t|)"
  ) +
  theme_minimal()



## ------------------------------------------------------------------
## 1  Overlay original mean curves … but strip all titles / axis labels
## ------------------------------------------------------------------

p_linear_comb <- p_linear +
  geom_line(data = df1, aes(x = x, y = y1),
            colour = "red",  linewidth = 1.4) +               # red solid
  geom_line(data = df1, aes(x = x, y = y2),
            colour = "blue", linetype = "dashed",
            linewidth = 1.4) +
  labs(title = NULL, x = NULL, y = NULL) +                    # remove titles
  theme_minimal() +
  theme(legend.position = "none",
        axis.title       = element_blank(),
        plot.title       = element_blank())

p_abs_comb <- p_abs +
  geom_line(data = df2, aes(x = x, y = y2),
            colour = "blue", linetype = "dashed",
            linewidth = 1.4) +                                # blue dashed
  geom_line(data = df2, aes(x = x, y = y1),
            colour = "red",  linewidth = 1.4) +               # red solid
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title       = element_blank(),
        plot.title       = element_blank())

## ------------------------------------------------------------------
## 2  Save at exactly 1026 × 673 px (dpi 100 avoids “squish”)
## ------------------------------------------------------------------

ggsave("linear_means_with_points.png",
       plot   = p_linear_comb,
       width  = 1026, height = 673,
       units  = "px", dpi = 100, bg = "white")

ggsave("abs_means_with_points.png",
       plot   = p_abs_comb,
       width  = 1026, height = 673,
       units  = "px", dpi = 100, bg = "white")


### Bigger labels: 

## ── tweak axis labels & font sizes ────────────────────────────────────────
big_axis_text  <- 18   # tick labels
big_axis_title <- 18   # “Time” label

p_linear_comb <- p_linear_comb +             # reuse existing object
  labs(x = "Time", y = NULL) +
  theme(
    axis.text  = element_text(size = big_axis_text),
    axis.title = element_text(size = big_axis_title, face = "plain")
  )

p_abs_comb <- p_abs_comb +
  labs(x = "Time", y = NULL) +
  theme(
    axis.text  = element_text(size = big_axis_text),
    axis.title = element_text(size = big_axis_title, face = "plain")
  )

## ── save (dimensions & dpi unchanged) ─────────────────────────────────────
ggsave("linear_means_with_points.png",
       plot   = p_linear_comb,
       width  = 1026, height = 673,
       units  = "px", dpi = 100, bg = "white")

ggsave("abs_means_with_points.png",
       plot   = p_abs_comb,
       width  = 1026, height = 673,
       units  = "px", dpi = 100, bg = "white")



