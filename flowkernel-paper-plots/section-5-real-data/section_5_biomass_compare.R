library(flowkernel)
library(mclust)
library(tidyverse)
library(ggplot2)
library(plotly)
library(readxl)

set.seed(123)
### PLOTTING FUNCTIONS ###

# Define a function to plot the gated biomass data
plot_manual_gated_biomass <- function(df) {
  # Define colors for each cluster
  line_colors <- c("Prochlorococcus" = "green", "Synechococcus" = "orange", "Picoeukaryotes" = "blue")
  
  # Create the ggplot with multiple line plots
  pi_plt <- ggplot(df, aes(x = time)) +
    geom_line(aes(y = Prochlorococcus, color = "Prochlorococcus")) +
    geom_line(aes(y = Synechococcus, color = "Synechococcus")) +
    geom_line(aes(y = Picoeukaryotes, color = "Picoeukaryotes")) +
    labs(x = "Time", y = "c_per_uL") +
    ggtitle("Biomass by Date") +
    scale_color_manual(name = "Cluster", values = line_colors)
  
  # Convert ggplot to plotly for interactivity
  pi_plotly <- ggplotly(pi_plt, dynamicTicks = TRUE)
  
  # Return the interactive plot
  return(pi_plotly)
}

plot_biomass_dates <- function(biomass, resp, dates) {
  K <- ncol(resp[[1]])  # Number of clusters
  ntimes <- length(resp)  # Number of time points
  
  # Check if the length of dates matches the number of time points
  if (length(dates) != ntimes) {
    stop("Length of 'dates' must match the number of time points in 'resp' and 'biomass'.")
  }
  
  # Initialize a list to hold data frames for each cluster
  data_list <- vector("list", K)
  
  for (k in 1:K) {
    # Compute biomass for each cluster at each time point
    cluster_biomass <- sapply(1:ntimes, function(tt) sum(resp[[tt]][, k] * biomass[[tt]]))
    
    # Create a data frame using the actual dates for the x-axis
    data_list[[k]] <- data.frame(time = dates, Cluster = paste("Cluster", k), Biomass = cluster_biomass)
  }
  
  # Combine all the data frames into one
  df <- do.call(rbind, data_list)
  
  # Create the ggplot with the actual dates on the x-axis
  plt <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = Biomass, color = Cluster)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Date", y = "Cluster Biomass") +
    ggplot2::ggtitle("Cluster Biomass Over Time") +
    ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d %H:%M", date_breaks = "1 day") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Convert ggplot to plotly for interactivity
  pi_plotly <- plotly::ggplotly(plt, dynamicTicks = TRUE)
  
  return(pi_plotly)
}


plot_combined_multiple_with_dates <- function(fit, biomass_data, cluster_nums, population_name, biomass_list, fit_dates, biomass_dates) {
  # Verify cluster numbers are provided
  if (length(cluster_nums) < 1) {
    stop("Please provide at least one cluster number.")
  }
  
  # Calculate the sum of responses across all specified clusters
  cluster_data <- sapply(seq_len(length(fit$resp)), function(tt) {
    sum(sapply(cluster_nums, function(cluster_num) {
      fit$resp[[tt]][, cluster_num] * biomass_list[[tt]]
    }))
  })
  
  # Ensure the cluster_data and fit_dates have the same length
  if (length(cluster_data) != length(fit_dates)) {
    stop("Mismatch between the length of cluster data and fit dates.")
  }
  
  # Align the actual biomass data based on population_name
  population_data <- biomass_data[[population_name]]
  
  # Create data frame for plotting
  df_combined <- data.frame(
    Date = c(fit_dates, biomass_dates),
    Biomass = c(cluster_data, population_data),
    Type = c(rep("Fitted Cluster Sum", length(cluster_data)), rep("Manual Gated Biomass", length(population_data)))
  )
  
  # Update ggplot title
  cluster_nums_str <- paste(cluster_nums, collapse = ", ")
  
  # Plot with ggplot
  plt <- ggplot(df_combined, aes(x = Date, y = Biomass, color = Type)) +
    geom_line(data = df_combined[df_combined$Type == "Fitted Cluster Sum", ], aes(x = Date, y = Biomass), linetype = "solid") +
    geom_line(data = df_combined[df_combined$Type == "Manual Gated Biomass", ], aes(x = Date, y = Biomass), linetype = "dashed") +
    labs(x = "Date", y = "Biomass") +
    ggtitle(paste("Comparison of Fitted Cluster Sum (Clusters:", cluster_nums_str, ") with Manual Gated Biomass for", population_name)) +
    scale_color_manual(values = c("Fitted Cluster Sum" = "red", "Manual Gated Biomass" = "blue"))
  
  # Convert to plotly for interactivity
  plotly_combined <- ggplotly(plt, dynamicTicks = TRUE)
  
  return(plotly_combined)
}

plot_combined_multiple <- function(fit,
                                   biomass_data,
                                   cluster_nums,
                                   population_name,
                                   biomass_list,
                                   fit_dates) {
  
  ## 1 ─ fitted-cluster series (numeric scalar per time-point)
  fitted_df <- tibble::tibble(
    Date    = as.POSIXct(fit_dates, tz = "UTC"),
    Biomass = purrr::map_dbl(seq_along(fit$resp), function(tt) {
      sum(fit$resp[[tt]][ , cluster_nums, drop = FALSE] *
            biomass_list[[tt]])
    }),
    Type    = "Fitted Cluster Sum"
  )
  
  ## 2 ─ manual series (works for wide or long formats)
  time_col <- intersect(names(biomass_data), c("date", "time", "Date", "Time"))[1]
  
  if (population_name %in% names(biomass_data)) {          # wide
    manual_df <- biomass_data |>
      dplyr::transmute(Date    = .data[[time_col]],
                       Biomass = as.numeric(.data[[population_name]]),
                       Type    = "Manual Gated Biomass")
  } else if (all(c("population", "c_per_uL") %in% names(biomass_data))) {  # long
    manual_df <- biomass_data |>
      dplyr::filter(population == population_name) |>
      dplyr::transmute(Date    = .data[[time_col]],
                       Biomass = as.numeric(c_per_uL),
                       Type    = "Manual Gated Biomass")
  } else {
    stop("`population_name` not found in manual data.")
  }
  
  ## 3 ─ plot
  combined_df <- dplyr::bind_rows(fitted_df, manual_df)
  
  ggplotly(
    ggplot2::ggplot(combined_df,
                    ggplot2::aes(Date, Biomass,
                                 colour = Type, linetype = Type)) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "Date",
                    y = "Biomass (c_per_uL)",
                    title = sprintf("Clusters %s vs manual %s",
                                    paste(cluster_nums, collapse = ", "),
                                    population_name)) +
      ggplot2::scale_colour_manual(
        values = c("Fitted Cluster Sum"   = "red",
                   "Manual Gated Biomass" = "blue")) +
      ggplot2::theme_minimal(),
    dynamicTicks = TRUE
  )
}


# Define a function to plot biomass over time
plot_biomass_over_time <- function(data) {
  # Define colors for each population
  line_colors <- c(
    "Allother" = "red",
    "Picoeukaryotes" = "blue",
    "Prochlorococcus" = "green",
    "Synechococcus" = "orange"
  )
  
  # Create the ggplot
  pi_plt <- ggplot(data, aes(x = date, y = c_per_uL, color = population)) +
    geom_line() +
    labs(x = "Date", y = "Biomass (c_per_uL)") +
    ggtitle("Biomass over Time for Different Populations") +
    scale_color_manual(name = "Population", values = line_colors) +
    theme_minimal()
  
  # Convert ggplot to plotly for interactivity
  pi_plotly <- ggplotly(pi_plt, dynamicTicks = TRUE)
  
  # Return the interactive plot
  return(pi_plotly)
}
### END PLOTTING FUNCTIONS DEFINITION


### Load manual data first
MGL1704_gated_data <- read.csv("~/real_data_analysis/MGL1704.csv")
# Parse the 'date' column to date-time objects
MGL1704_gated_data$date <- as.POSIXct(MGL1704_gated_data$date, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
plot_biomass_over_time(MGL1704_gated_data)


### Load ungated data to run the algorithm on ###
grid_file_path <- "~/SeaFlow_Data/2022-11-01-diam-chl-pe-Qc-PSD-combined-popcycle-4.15.0-withoutBoundaryPoints.psd-grid.parquet"

# Read the parquet file
grid <- arrow::read_parquet(grid_file_path)

data_file_path <- "~/SeaFlow_Data/2022-11-01-diam-chl-pe-Qc-PSD-combined-popcycle-4.15.0-withoutBoundaryPoints.psd-hourly.parquet"
psd_hourly <- arrow::read_parquet(data_file_path) 

psd_hourly <- psd_hourly %>%
  mutate(diam = grid$diam[psd_hourly$diam_coord],
         Qc = grid$Qc[psd_hourly$Qc_coord],
         chl_small = grid$chl_small[psd_hourly$chl_small_coord],
         pe = grid$pe[psd_hourly$pe_coord]) %>%
  select(-c(diam_coord, Qc_coord, chl_small_coord, pe_coord))

MGL1704_data <- psd_hourly[psd_hourly$cruise == "MGL1704", ]

MGL1704_y_list <- MGL1704_data %>%
  group_by(date) %>%
  group_map(~ .x %>% select(diam, chl_small, pe)) %>%
  map(as.matrix)

dates_list <- MGL1704_data %>%
  group_by(date) %>%
  group_keys() %>%
  pull(date)

MGL1704_y_list <- list(date = dates_list, data = MGL1704_y_list)

MGL1704_biomass_list <- MGL1704_data %>% 
  group_by(date) %>%
  group_map(~ .x) %>% 
  map(~ .x %>% select(Qc_sum_per_uL))

# Convert each data frame to a numeric vector
MGL1704_biomass_list <- map(MGL1704_biomass_list, function(df) {
  # Extract the column from the data frame and convert it to numeric
  numeric_vector <- as.numeric(unlist(df))
  return(numeric_vector)
})

#Try log scale

MGL1704_y_list_log <- lapply(MGL1704_y_list$data, function(mat) log(mat))
set.seed(123)
MGL1704_fit_7_log <- kernel_em(MGL1704_y_list_log, K = 7, hmu = 23, hSigma = 15, hpi = 23, dates = MGL1704_y_list$date,
                               biomass = MGL1704_biomass_list)
plot_data(MGL1704_y_list_log)
plot_data_and_model(y = MGL1704_y_list_log, z =  MGL1704_fit_7_log$zest, mu = MGL1704_fit_7_log$mu)

set.seed(123)
MGL1704_fit_8_log <- kernel_em(MGL1704_y_list_log, K = 8, hmu = 23, hSigma = 15, hpi = 23, dates = MGL1704_y_list$date,
                               biomass = MGL1704_biomass_list)

plot_combined_multiple <- function(fit, biomass_data, cluster_nums, population_name,
                                   biomass_list, fit_dates,
                                   interactive = FALSE, save = FALSE) {
  
  ## ── fitted series ─────────────────────────────────────────────
  fitted_df <- tibble::tibble(
    Date    = as.POSIXct(fit_dates, tz = "UTC"),
    Biomass = purrr::map_dbl(seq_along(fit$resp), \(tt)
                             sum(fit$resp[[tt]][ , cluster_nums, drop = FALSE] *
                                   biomass_list[[tt]])),
    Line    = "Fitted"
  )
  
  ## ── manual series (wide or long) ─────────────────────────────
  time_col <- intersect(names(biomass_data),
                        c("date", "time", "Date", "Time"))[1]
  
  if (population_name %in% names(biomass_data)) {          # wide
    manual_df <- biomass_data |>
      dplyr::transmute(Date    = .data[[time_col]],
                       Biomass = as.numeric(.data[[population_name]]),
                       Line    = "Manual")
  } else {                                                 # long
    manual_df <- biomass_data |>
      dplyr::filter(population == population_name) |>
      dplyr::transmute(Date    = .data[[time_col]],
                       Biomass = as.numeric(c_per_uL),
                       Line    = "Manual")
  }
  
  df <- dplyr::bind_rows(fitted_df, manual_df)
  
  ## ── plot ─────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(Date, Biomass,
                                        colour = Line, linetype = Line)) +
    ggplot2::geom_line(size = 0.6) +
    ggplot2::scale_colour_manual(values = c(Fitted = "red", Manual = "blue"),
                                 guide = "none") +
    ggplot2::scale_linetype_manual(values = c(Fitted = "solid", Manual = "dashed"),
                                   guide = "none") +
    ggplot2::scale_x_datetime(breaks = scales::breaks_pretty(5)) +
    ggplot2::labs(x = NULL, y = "Biomass (c per µL)") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      panel.border       = ggplot2::element_rect(colour = "black",
                                                 fill = NA, linewidth = 0.3),
      axis.text          = ggplot2::element_text(size = 10),
      panel.grid.major   = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank()
    )
  
  ## ── save, if requested ───────────────────────────────────────
  if (save) {
    pop_clean  <- stringr::str_replace_all(tolower(population_name), "\\s+", "-")
    clust_part <- paste(cluster_nums, collapse = "-")
    fname      <- sprintf("biomass-compare-%s-%s.png", pop_clean, clust_part)
    ggplot2::ggsave(fname, p, width = 6, height = 3.5, dpi = 300)
  }
  
  ## ── return object ────────────────────────────────────────────
  if (interactive) {
    return(plotly::ggplotly(p, dynamicTicks = TRUE))
  }
  p
}



# Compare manual and kernel-em
plot_biomass_dates(biomass = MGL1704_biomass_list, resp = MGL1704_fit_8_log$resp, dates =MGL1704_y_list$date)
plot_biomass_over_time(MGL1704_gated_data)
plot_manual_gated_biomass(MGL1704_data_gated_subset)
plot_manual_gated_biomass(MGL1704_gated_data)
plot_biomass_over_time(MGL1704_gated_data)

prochloro_compare <- plot_combined_multiple(fit         = MGL1704_fit_8_log,
                                  biomass_data = MGL1704_gated_data,
                                  cluster_nums = c(1, 8),
                                  population_name = "Prochlorococcus",
                                  biomass_list  = MGL1704_biomass_list,
                                  fit_dates     = MGL1704_y_list$date,
                                  save = TRUE)

synecho_compare <- plot_combined_multiple(fit         = MGL1704_fit_8_log,
                       biomass_data = MGL1704_gated_data,
                       cluster_nums = c(5),
                       population_name = "Synechococcus",
                       biomass_list  = MGL1704_biomass_list,
                       fit_dates     = MGL1704_y_list$date,
                       save = TRUE)
pico_compare <- plot_combined_multiple(fit         = MGL1704_fit_8_log,
                       biomass_data = MGL1704_gated_data,
                       cluster_nums = c(6),
                       population_name = "Picoeukaryotes",
                       biomass_list  = MGL1704_biomass_list,
                       fit_dates     = MGL1704_y_list$date,
                       save = TRUE)
all_other_compare <- plot_combined_multiple(fit         = MGL1704_fit_8_log,
                       biomass_data = MGL1704_gated_data,
                       cluster_nums = c(2, 3, 4, 7),
                       population_name = "Allother",
                       biomass_list  = MGL1704_biomass_list,
                       fit_dates     = MGL1704_y_list$date,
                       save = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# tweak_and_save() ------------------------------------------------------------
# p            : ggplot from plot_combined_multiple()
# filename     : file name for ggsave()
# add_legend   : TRUE  → keep a legend (“Manual gated”, “Kernel-EM gated”)
#                FALSE → no legend (default, mimics your current plots)
# width, height: physical size in inches; dpi: as usual
# ─────────────────────────────────────────────────────────────────────────────
tweak_and_save <- function(p,
                           filename,
                           add_legend = FALSE,
                           width  = 6,
                           height = 3.5,
                           dpi    = 300) {
  
  ## ── base clean-ups ───────────────────────────────────────────────
  p2 <- p +
    theme_classic(base_size = 12) +                       # keep white bg
    theme(
      # axis titles a bit larger
      axis.title = element_text(size = 13),
      # light grey grid on white
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
      # single uniform border; remove redundant axis lines
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.3),
      axis.line        = element_blank(),
      # legend handled below
      legend.position  = "none"
    )
  
  ## ── optional legend ─────────────────────────────────────────────
  if (add_legend) {
    p2 <- p2 +
      ## correct colour / linetype and put Kernel-EM first -----------
    scale_colour_manual(
      breaks = c("Fitted", "Manual"),                     # order in legend
      values = c(Fitted = "red", Manual = "blue"),
      labels = c("Kernel-EM gated", "Manual gated"),
      name   = NULL
    ) +
      scale_linetype_manual(
        breaks = c("Fitted", "Manual"),
        values = c(Fitted = "solid", Manual = "dashed"),
        labels = c("Kernel-EM gated", "Manual gated"),
        name   = NULL
      ) +
      guides(colour = guide_legend(override.aes = list(size = 0.9))) +
      theme(
        legend.position      = c(1, 1),   # inside top-right
        legend.justification = c(1, 1),
        legend.background    = element_rect(fill = alpha("white", 0.8),
                                            colour = NA),
        legend.key.width     = unit(1.1, "lines")
      )
  }
  
  
  ## ── save ─────────────────────────────────────────────────────────
  ggsave(filename, p2, width = width, height = height, dpi = dpi)
  
  invisible(p2)     # return the tweaked plot invisibly
}

# ─────────────────────────────────────────────────────────────────
# EXAMPLE USAGE
# tweak without legend
tweak_and_save(all_other_compare, "all_other_vs_rest.png")

# tweak WITH legend (inside the frame)
tweak_and_save(synecho_compare, "synec_vs_cluster5_withlegend.png",
               add_legend = TRUE)


