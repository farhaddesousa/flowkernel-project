library(ggplot2)
#load real data
library(arrow)
library(tidyverse)
library(magrittr)
library(dplyr)
library(flowkernel)

#### PARTICLE LEVEL DATA LOAD AND FIT ### 
#particle gated gated MGL for pico comparison
## Try all data
library(arrow)
# Define the zip file path and output directory
zipfile <- "MGL1704_vct.zip"
output_dir <- "extracted_MGL1704"

# List the contents of the zip file
zip_content <- unzip(zipfile, list = TRUE)

# Filter the .parquet files in the MGL1704 directory
parquet_files <- grep("MGL1704_vct/.*\\.parquet$", zip_content$Name, value = TRUE)

# Unzip only the .parquet files into the output directory
unzip(zipfile, files = parquet_files, exdir = output_dir)


# List all the .parquet files in the extracted directory
parquet_files <- list.files(
  "extracted_MGL1704/MGL1704_vct",
  pattern = "\\.parquet$",
  full.names = TRUE
)

# Read the .parquet files into a list of data frames
MGL1704_particle_gated <- lapply(parquet_files, read_parquet)

# View the first data frame
View(MGL1704_particle_gated[[1]])
class(MGL1704_particle_gated[[1]])

# Extract only the "fsc_small", "chl_small", and "pe" columns
particle_data <- lapply(MGL1704_particle_gated, function(df) df[, c("fsc_small", "chl_small", "pe")])
particle_data <- lapply(particle_data, function(df) as.matrix(df))


View(particle_data[[1]])
length(particle_data)


#Try from package
particle_fit <- flowkernel::kernel_em(y = particle_data, K = 10, hmu = 20, hSigma = 15, hpi = 20)

#### BINNED DATA LOADING ###
grid_file_path <- "~/SeaFlow_Data/2022-11-01-diam-chl-pe-Qc-PSD-combined-popcycle-4.15.0-withoutBoundaryPoints.psd-grid.parquet"

# Read the parquet file
grid <- arrow::read_parquet(grid_file_path)

data_file_path <- "~/SeaFlow_Data/2022-11-01-diam-chl-pe-Qc-PSD-combined-popcycle-4.15.0-withoutBoundaryPoints.psd-hourly.parquet"
options(arrow.unsafe_metadata = TRUE)
psd_hourly <- arrow::read_parquet(data_file_path)


psd_hourly <- psd_hourly %>%
  mutate(diam = grid$diam[psd_hourly$diam_coord],
         Qc = grid$Qc[psd_hourly$Qc_coord],
         chl_small = grid$chl_small[psd_hourly$chl_small_coord],
         pe = grid$pe[psd_hourly$pe_coord]) %>%
  select(-c(diam_coord, Qc_coord, chl_small_coord, pe_coord))

MGL1704_data <- psd_hourly[psd_hourly$cruise == "MGL1704", ]


### NEW FORMAT MAR-3-25
# Subset data for a specific cruise
MGL1704_data <- psd_hourly[psd_hourly$cruise == "MGL1704", ]

# Group measurements by date and convert to matrix
MGL1704_y_list <- MGL1704_data %>%
  group_by(date) %>%
  group_map(~ .x %>% select(diam, chl_small, pe) %>% as.matrix)

# Extract unique dates
dates_list <- MGL1704_data %>%
  group_by(date) %>%
  group_keys() %>%
  pull(date)

# Gather biomass measurements by date
MGL1704_biomass_list <- MGL1704_data %>% 
  group_by(date) %>%
  group_map(~ .x %>% select(Qc_sum_per_uL)) %>%
  map(~ as.numeric(unlist(.x)))

# Log-transform data (optional, depending on your analysis)
MGL1704_log_y_list <- lapply(MGL1704_y_list, function(mat) {
  if(any(mat <= 0)) {
    warning("Some values are non-positive. log() may yield -Inf or NaN.")
  }
  log(mat)
})

# Combine processed data into a single list
MGL1704_lists <- list(
  y     = MGL1704_y_list,
  dates    = dates_list,
  biomass  = MGL1704_biomass_list,
  log_y = MGL1704_log_y_list
)

dates_trim <- head(dates_list, dim(particle_fit$mu)[1])   # keep first 355 dates


## PLOT MEANS ##
library(ggplot2)
library(dplyr)
library(tidyr)

# ══════════════════════════════════════════════════════════════════════
# Publish-ready plot of selected cluster means (1 dimension)
# ══════════════════════════════════════════════════════════════════════
plot_cluster_means <- function(mu,
                               dates_vec,
                               clusters   = c(3, 5, 8),
                               dim        = 1,
                               y          = NULL,
                               line_width = 0.5) {          # ← new arg
  
  ## ---- basic checks --------------------------------------------------
  if (length(dim(mu)) != 3L)
    stop("'mu' must be a 3-D array (T × K × d)")
  if (!inherits(dates_vec, "POSIXt"))
    dates_vec <- as.POSIXct(dates_vec, tz = "UTC")   # convert if user gave strings
  
  T  <- dim(mu)[1]          # number of time points
  d  <- dim(mu)[3]
  
  if (dim < 1 || dim > d)
    stop("'dim' out of range (must be 1 … ", d, ")")
  
  ## ---- align lengths -------------------------------------------------
  if (length(dates_vec) != T) {
    warning("dates_vec length (", length(dates_vec),
            ") ≠ number of time slices in mu (", T,
            "); using the first T timestamps.")
    dates_vec <- dates_vec[seq_len(T)]
  }
  
  ## ---- build tidy tibble --------------------------------------------
  slice        <- mu[ , clusters, dim]                  # T × |clusters|
  colnames(slice) <- paste0("cl", clusters)
  
  df <- as.data.frame(slice) |>
    mutate(Time = dates_vec) |>
    pivot_longer(-Time, names_to = "Cluster", values_to = "Mean")
  
  ## ---- y-axis label --------------------------------------------------
  ylab <- if (!is.null(y) && !is.null(colnames(y[[1]]))) {
    colnames(y[[1]])[dim]
  } else {
    paste0("Dimension ", dim)
  }
  
  ## ---- colours -------------------------------------------------------
  cols <- scales::hue_pal()(length(clusters))
  
  ## ---- ggplot --------------------------------------------------------
  p <- ggplot(df, aes(Time, Mean,
                      colour = Cluster, group = Cluster)) +
    geom_line(linewidth = line_width) +        # ← use it here
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_datetime(date_labels = "%b %d", date_breaks = "4 days") +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.margin      = margin(5, 10, 5, 5)
    )
  
  p
}

# ── example call ───────────────────────────────────────────────────────
dates_trim <- head(dates_list, dim(particle_fit$mu)[1])   # if you still need trimming
plot_cluster_means(particle_fit$mu,
                   dates_vec = dates_trim,
                   clusters  = c(3, 5, 8),
                   dim       = 1,
                   y         = particle_data)

library(ggplot2)

# ── tweak axis text and save ──────────────────────────────────────────
save_cluster_means <- function(p,
                               file,
                               width  = 4,     # inches — three side-by-side ≈ full textwidth
                               height = 2.4,
                               dpi    = 300) {
  
  p_out <- p +
    theme(
      axis.title = element_text(size = 12),   # larger labels
      axis.text  = element_text(size = 11),                  # slightly larger ticks
      legend.position = "none"                               # keep it clean
    )
  
  ggsave(filename = file,
         plot     = p_out,
         width    = width,
         height   = height,
         units    = "in",
         dpi      = dpi,
         bg       = "white")
  
  invisible(p_out)   # return tweaked plot (useful in notebooks)
}

# ── example usage ────────────────────────────────────────────────────
p_raw <- plot_cluster_means(particle_fit$mu,
                            dates_vec = dates_trim,
                            clusters  = c(3, 5, 8),
                            dim       = 1,
                            y         = particle_data,
                            line_width = 0.5)     # thinner lines

save_cluster_means(p_raw, "pico_cluster_means_dim1.png")



## END PLOT MEANS ## 

## CELL ABUNDANCE BY RESPONSIBILITIES ##
# ─────────────────────────────────────────────────────────────────────────────
# plot_cluster_abundance() ----------------------------------------------------
#  • resp_list  : particle_fit$resp  (list of length T, each item n_i × K)
#  • dates_vec  : same POSIXct vector you used for means
#  • clusters   : which columns (clusters) to plot
#  • save_path  : file name – if NULL the plot is returned but not written
# ─────────────────────────────────────────────────────────────────────────────
plot_cluster_abundance <- function(resp_list,
                                   dates_vec,
                                   clusters  = c(3, 5, 8),
                                   line_width = 0.5,     # <── add this
                                   save_path  = NULL,
                                   width  = 4,
                                   height = 2.4,
                                   dpi    = 300) {
  
  ## ---- consistency checks -----------------------------------------------
  T <- length(resp_list)
  if (length(dates_vec) != T)
    stop("dates_vec length (", length(dates_vec),
         ") ≠ length(resp_list) (", T, ").")
  
  K <- ncol(resp_list[[1]])
  if (any(clusters < 1) || any(clusters > K))
    stop("clusters must be within 1 … ", K)
  
  ## ---- soft counts per cluster & time ------------------------------------
  # matrix: rows = time points, cols = selected clusters
  counts_mat <- sapply(clusters, function(k)
    vapply(resp_list, function(mat) sum(mat[ , k]), numeric(1)))
  
  colnames(counts_mat) <- paste0("cl", clusters)
  
  ## ---- long data frame ---------------------------------------------------
  df <- as.data.frame(counts_mat) |>
    dplyr::mutate(Time = as.POSIXct(dates_vec, tz = "UTC")) |>
    tidyr::pivot_longer(-Time, names_to = "Cluster", values_to = "Abundance")
  
  ## ---- colours -----------------------------------------------------------
  cols <- scales::hue_pal()(length(clusters))
  
  ## ---- ggplot ------------------------------------------------------------
  p <- ggplot(df, aes(Time, Abundance,
                      colour = Cluster, group = Cluster)) +
    geom_line(linewidth = line_width) +        # <── use it here
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_datetime(date_labels = "%b %d", date_breaks = "4 days") +
    labs(x = NULL, y = "Cell Abundance") +
    theme_minimal(base_size = 13) +
    theme(
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      plot.margin      = margin(5, 10, 5, 5)
    )
  
  if (!is.null(save_path))
    ggsave(save_path, p, width = width, height = height,
           units = "in", dpi = dpi, bg = "white")
  
  invisible(p)
}

# ── example ---------------------------------------------------------------
abund_plot <- plot_cluster_abundance(
  resp_list  = particle_fit$resp,
  dates_vec  = dates_trim,
  clusters   = c(3, 5, 8),
  save_path  = "pico_cluster_abundance.png"     # omit or set NULL to skip saving
)

library(scales)  # for nicer tick labels (optional)

abund_plot_log <- abund_plot +
  scale_y_log10(labels = comma_format()) +   # log scale, ticks still in raw units
  labs(y = "Cell Abundance")                 # keep original axis title

ggsave("pico_abundance_log.png",
       abund_plot_log, width = 4, height = 2.4,
       units = "in", dpi = 300, bg = "white")



## END CELL ABUNDANCE RESPONSIBILITIES##

## CELL ABUNDANCE BY PI ##
##  π-weighted CELL ABUNDANCE  ──────────────────────────────────────────────
plot_pi_abundance <- function(pi_mat,             # T × K matrix
                              data_list,          # list of length T, each item n_i × d
                              dates_vec,          # POSIXct vector of length T
                              clusters   = c(3, 5, 8),
                              line_width = 0.5,
                              save_path  = NULL,
                              width  = 4,
                              height = 2.4,
                              dpi    = 300) {
  
  ## ---- consistency checks ------------------------------------------------
  T <- nrow(pi_mat)
  if (length(dates_vec) != T)
    stop("dates_vec length (", length(dates_vec), ") ≠ nrow(pi_mat) (", T, ").")
  
  K <- ncol(pi_mat)
  if (any(clusters < 1) || any(clusters > K))
    stop("clusters must be within 1 … ", K)
  
  if (length(data_list) != T)
    stop("data_list length must equal nrow(pi_mat).")
  
  ## ---- π-weighted counts per cluster & time ------------------------------
  n_vec <- vapply(data_list, nrow, integer(1))            # total cells each time
  counts_mat <- sapply(clusters, function(k) pi_mat[, k] * n_vec)
  colnames(counts_mat) <- paste0("cl", clusters)
  
  ## ---- long data frame ---------------------------------------------------
  df <- as.data.frame(counts_mat) |>
    dplyr::mutate(Time = as.POSIXct(dates_vec, tz = "UTC")) |>
    tidyr::pivot_longer(-Time, names_to = "Cluster", values_to = "Abundance")
  
  ## ---- colours -----------------------------------------------------------
  cols <- scales::hue_pal()(length(clusters))
  
  ## ---- ggplot ------------------------------------------------------------
  p <- ggplot(df, aes(Time, Abundance, colour = Cluster, group = Cluster)) +
    geom_line(linewidth = line_width) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_datetime(date_labels = "%b %d", date_breaks = "4 days") +
    labs(x = NULL, y = "Cell Abundance") +
    theme_minimal(base_size = 13) +
    theme(
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      plot.margin      = margin(5, 10, 5, 5)
    )
  
  if (!is.null(save_path))
    ggsave(save_path, p, width = width, height = height,
           units = "in", dpi = dpi, bg = "white")
  
  invisible(p)
}

pi_abund_plot <- plot_pi_abundance(
  pi_mat    = particle_fit$pi,
  data_list = particle_data,
  dates_vec = dates_trim,
  clusters  = c(3, 5, 8),
  save_path = "pico_cluster_abundance_pi.png"
)

library(scales)
pi_abund_plot_log <- pi_abund_plot +
  scale_y_log10(labels = comma_format()) +
  labs(y = "Cell Abundance")

ggsave("pico_abundance_pi_log.png",
       pi_abund_plot_log, width = 4, height = 2.4,
       units = "in", dpi = 300, bg = "white")



## END CELL ABUNDANCE BY PI ## 

View(particle_fit$resp[[1]])

################## Next Cruise KM1713 ###########################
#prepare gated data
KM1713_data_gated <- gated_data[gated_data$cruise == "KM1713", ]
KM1713_data_gated <- KM1713_data_gated %>%
  # Round down time to the nearest hour
  mutate(time = floor_date(time, unit = "hour")) %>%
  # Group by the new hourly time and cruise (if necessary)
  group_by(time) %>%
  # Calculate the mean of all numeric columns
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop')
head(KM1713_data_gated)
common_times_KM1713 <- intersect(KM1713_data_gated$time, KM1713_y_list$date)

# Subset the columns of interest
KM1713_data_gated_subset <- KM1713_data_gated[, c("time", "biomass_prochloro", "biomass_synecho", "biomass_picoeuk")]
View(KM1713_data_gated_subset)
# Rename the columns for better readability (optional)
colnames(KM1713_data_gated_subset) <- c("time", "Prochlorococcus", "Synechococcus", "Picoeukaryotes")

# Call the function to create and display the plot
plot_manual_gated_biomass(KM1713_data_gated_subset)

#Prepare data for algo:
KM1713_data <- psd_hourly[psd_hourly$cruise == "KM1713", ]

KM1713_y_list <- KM1713_data %>%
  group_by(date) %>%
  group_map(~ .x %>% select(diam, chl_small, pe)) %>%
  map(as.matrix)

dates_list <- KM1713_data %>%
  group_by(date) %>%
  group_keys() %>%
  pull(date)

KM1713_y_list <- list(date = dates_list, data = KM1713_y_list)

KM1713_biomass_list <- KM1713_data %>% 
  group_by(date) %>%
  group_map(~ .x) %>% 
  map(~ .x %>% select(Qc_sum_per_uL))

# Convert each data frame to a numeric vector
KM1713_biomass_list <- map(KM1713_biomass_list, function(df) {
  # Extract the column from the data frame and convert it to numeric
  numeric_vector <- as.numeric(unlist(df))
  return(numeric_vector)
})

#Try log scale

KM1713_y_list_log <- lapply(KM1713_y_list$data, function(mat) log(mat))

KM1713_fit_7_log <- kernel_em_dates(KM1713_y_list_log, K = 7, hmu = 23, hSigma = 15, hpi = 23, dates = KM1713_y_list$date,
                                    biomass = KM1713_biomass_list)
KM1713_fit_10_log <- kernel_em_dates(KM1713_y_list_log, K = 10, hmu = 23, hSigma = 15, hpi = 23, dates = KM1713_y_list$date,
                                     biomass = KM1713_biomass_list)
# Compare manual and kernel-em K = 7
plot_biomass_dates(biomass = KM1713_biomass_list, resp = KM1713_fit_7_log$resp, dates =KM1713_y_list$date)
plot_manual_gated_biomass(KM1713_data_gated_subset)

plot_combined_multiple_with_dates(
  fit = KM1713_fit_7_log, 
  biomass_data = KM1713_data_gated_subset, 
  cluster_nums = c(1),  
  population_name = "Picoeukaryotes",
  biomass_list = KM1713_biomass_list,
  fit_dates = KM1713_y_list$date,
  biomass_dates = KM1713_data_gated_subset$time
)

# Compare manual and kernel-em K = 10
plot_biomass_dates(biomass = KM1713_biomass_list, resp = KM1713_fit_10_log$resp, dates =KM1713_y_list$date)
plot_manual_gated_biomass(KM1713_data_gated_subset)

plot_combined_multiple_with_dates(
  fit = KM1713_fit_10_log, 
  biomass_data = KM1713_data_gated_subset, 
  cluster_nums = c(1),  
  population_name = "Picoeukaryotes",
  biomass_list = KM1713_biomass_list,
  fit_dates = KM1713_y_list$date,
  biomass_dates = KM1713_data_gated_subset$time
)