#load real data
library(arrow)
library(tidyverse)
library(magrittr)
library(dplyr)
library(flowkernel)
library(mclust)
library(plo)
library(dplyr)
library(purrr)

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

MGL1704_y_list_log <- readRDS("MGL1704_y_list_log.rds")
MGL1704_fit_7_log  <- readRDS("MGL1704_fit_7_log.rds")

# raw_fig   <- flowkernel::plot_data(MGL1704_y_list_log)
# gated_fig <- flowkernel::plot_data_and_model(
#   y  = MGL1704_y_list_log,
#   z  = MGL1704_fit_7_log$zest,
#   mu = MGL1704_fit_7_log$mu
# )


library(htmlwidgets)
library(webshot2)

save_plotly_png <- function(fig, file, width, height, scale = 3){
  tmp <- tempfile(fileext = ".html")
  saveWidget(as_widget(fig), tmp, selfcontained = TRUE)
  webshot(tmp, file = file, vwidth = width, vheight = height, zoom = scale)
}

`%||%` <- function(a,b) if (is.null(a)) b else a

freeze_frame_hd <- function(p,
                            t_value,
                            file,
                            data_pt_size   = 5,
                            center_pt_size = 12,
                            axis_title_sz  = 36,
                            axis_tick_sz   = 24,
                            width = 1600, height = 1400, scale = 3,
                            cam_eye   = list(x = -1.05, y = 1.58, z = 0.45),
                            cam_up    = list(x = 0, y = 0, z = 1),
                            cam_center= list(x = 0, y = 0, z = 0),
                            # how far to push titles beyond the axis ranges (fraction of range)
                            offset_frac = 0.10,
                            # outer margins around the scene (px)
                            margins = list(l = 60, r = 40, b = 60, t = 25)) {
  
  `%||%` <- function(a,b) if (is.null(a)) b else a
  
  pb <- plotly::plotly_build(p)
  
  # choose frame by name; fallback to index
  fr_idx <- which(vapply(pb$x$frames, function(fr) fr$name, "") == as.character(t_value))
  if (!length(fr_idx)) fr_idx <- as.integer(t_value)
  fr <- pb$x$frames[[fr_idx]]
  
  # freeze the frame
  pb$x$data   <- fr$data
  pb$x$frames <- NULL
  pb$x$layout$title$text <- NULL
  pb$x$layout$showlegend  <- FALSE
  pb$x$layout$updatemenus <- NULL
  pb$x$layout$sliders     <- NULL
  
  # scene + camera
  sc <- pb$x$layout$scene %||% list()
  sc$camera <- list(eye = cam_eye, up = cam_up, center = cam_center)
  
  # keep tick sizes; turn off built-in titles (we'll add annotations)
  for (ax in c("xaxis","yaxis","zaxis")) {
    sc[[ax]] <- modifyList(sc[[ax]] %||% list(), list(
      title    = list(text = ""),
      tickfont = list(size = axis_tick_sz),
      ticklen  = 8,
      automargin = TRUE
    ))
  }
  
  # compute data ranges for precise placement
  get_rng <- function(coord) {
    vals <- unlist(lapply(fr$data, function(tr) if (identical(tr$type, "scatter3d")) tr[[coord]]))
    range(vals, na.rm = TRUE)
  }
  xr <- get_rng("x"); yr <- get_rng("y"); zr <- get_rng("z")
  dx <- offset_frac * diff(xr)
  dy <- offset_frac * diff(yr)
  dz <- offset_frac * diff(zr)
  
  # axis titles as 3D annotations placed just outside the cube
  sc$annotations <- list(
    # X title near positive-x corner, pushed outward
    list(x = xr[2] + dx, y = yr[1] - dy, z = zr[1] - 0.25*dz,
         xref="x", yref="y", zref="z", showarrow=FALSE, text="Diameter",
         font=list(size=axis_title_sz), xanchor="left"),
    # Y title near positive-y edge
    list(x = xr[1] - 0.25*dx, y = yr[2] + dy, z = zr[1] - 0.25*dz,
         xref="x", yref="y", zref="z", showarrow=FALSE, text="Chl",
         font=list(size=axis_title_sz), xanchor="center"),
    # Z title near positive-z edge
    list(x = xr[1] - 0.25*dx, y = yr[1] - 0.25*dy, z = zr[2] + dz,
         xref="x", yref="y", zref="z", showarrow=FALSE, text="Pe",
         font=list(size=axis_title_sz), xanchor="right")
  )
  
  sc$aspectmode  <- "manual"
  sc$aspectratio <- list(x = 1, y = 1, z = 1)
  pb$x$layout$scene  <- sc
  
  # tighter margins → bigger plot area in the PNG
  pb$x$layout$margin <- margins
  
  # readable markers at high DPI
  pb$x$data <- lapply(pb$x$data, function(tr){
    if (identical(tr$type, "scatter3d")) {
      tr$marker <- tr$marker %||% list()
      tr$marker$size    <- if (length(tr$x)==1L && length(tr$y)==1L && length(tr$z)==1L) center_pt_size else data_pt_size
      tr$marker$opacity <- tr$marker$opacity %||% 0.9
    }
    tr
  })
  
  # save (uses your save_plotly_png if defined; else webshot2)
  if (exists("save_plotly_png")) {
    save_plotly_png(pb, file, width = width, height = height, scale = scale)
  } else {
    tmp <- tempfile(fileext = ".html")
    htmlwidgets::saveWidget(plotly::as_widget(pb), tmp, selfcontained = TRUE)
    webshot2::webshot(tmp, file = file, vwidth = width, vheight = height, zoom = scale)
  }
  
  invisible(normalizePath(file))
}




`%||%` <- function(a,b) if (is.null(a)) b else a


# -------------------- export the six HD panels -------------------------
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

for (tt in c(1, 82, 346)) {
  freeze_frame_hd(raw_fig,   t_value = tt,
                  file = sprintf("Figures/MGLraw-t%d.png", tt))
  freeze_frame_hd(gated_fig, t_value = tt,
                  file = sprintf("Figures/MGL-gated-t%d.png", tt))
}


