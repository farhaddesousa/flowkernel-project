# This is the main simulator file
library(mclust)
library(ggplot2)
library(simulator) # this file was created under simulator version 0.2.4
library(flowkernel)
source("model_functions.R")
source("method_functions.R")
source("eval_functions.R")

# Disappearance simulation
set.seed(123)
hmu <- 20
hpi <- 10

disappearance_durations <- seq(0, 60, by = 5)
transition_duration <- 30
hSigma <- 100
n_runs <- 100
ncores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))

# Create the disappearance simulation object
sim <- new_simulation(name = "disappearance_study",
                      label = "Disappearance Study Simulation") %>%
  generate_model(make_disappearance_model,
                 transition_duration = transition_duration,
                 disappearance_duration = as.list(disappearance_durations),
                 vary_along = "disappearance_duration") %>%
  simulate_from_model(nsim = n_runs)
## 2) individual RDS files (handy if only one object is needed later)
saveRDS(sim,        file = "rdata/dis_sim_seed_123.rds",        compress = "xz")


method_sim <- sim %>% run_method(
  list(kernel_em_CV_method, const_method, hungarian_method),
  parallel = list(
    socket_names = ncores,    # Number of cores to use
    libraries = c("flowkernel", "foreach", "doParallel", "reshape2", "mclust")   # Libraries to load on each worker
  )
)

saveRDS(method_sim, file = "rdata/dis_method_sim_seed_123", compress = "xz")

full_sim <- method_sim %>% evaluate(list(combined_concatenated_rand_index_metric, combined_concatenated_ARI_metric))
saveRDS(full_sim,   file = "rdata/dis_full_sim_seed_123.rds",   compress = "xz")

# Save
## ─────────────────────────────────────────────────────────────────
## after full_sim is created … save everything you’ll need later
## ─────────────────────────────────────────────────────────────────

## 1) one all-in-one workspace file (easy to reload)
save(sim, method_sim, full_sim,
     file = "rdata/disappearance_workspace_seed_123.RData",
     compress = "xz")          # smallest size, still fast





## LOAD SAVED OBJECTS
## load packages FIRST (objects refer to their S3/S4 classes)

## Option A: pull everything back at once
load("rdata/disappearance_workspace.RData")
# objects sim, method_sim, full_sim now exist exactly as before

## Option B: grab just what you need
full_sim <- readRDS("rdata/full_sim.rds")


plot_rand_dis <- plot_eval_by(full_sim, "combined_concatenated_rand_index", varying = "disappearance_duration", main = "Cluster Disappearance Simulation", 
             xlab = "Disappearance Duration", ylab = "Rand Index")

plot_adj_rand_dis <- plot_eval_by(full_sim, "combined_concatenated_adjusted_rand_index", varying = "disappearance_duration", main = "Cluster Disappearance Simulation", 
                              xlab = "Disappearance Duration", ylab = "Adjusted Rand Index")

#Save paper formated plots:


#─────────────────────────────────────────────────────────────────────────────
# style_eval_plot()  –– tidy up a plot_eval_by() ggplot and (optionally) save
#─────────────────────────────────────────────────────────────────────────────
# p        : ggplot returned by plot_eval_by()
# file     : if non-NULL, write PNG here
# width/height/dpi : physical size for ggsave()
#─────────────────────────────────────────────────────────────────────────────         # for unit()
library(ggplot2)
style_eval_plot <- function(p,
                            show_legend = TRUE,
                            file        = NULL,
                            width  = 4,
                            height = 2.4,
                            dpi    = 300) {
  
  base_theme <- theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor = element_line(colour = "grey92", linewidth = 0.15),
      axis.title       = element_text(size = 10, face = "plain")
    )
  
  legend_theme <- if (show_legend) {
    theme(
      legend.position      = c(0, 0),                 # bottom-left inside
      legend.justification = c(0, 0),
      legend.background    = element_rect(fill = alpha("white", 0.75),
                                          colour = NA),
      legend.key.size      = unit(0.3, "lines"),      # smaller box
      legend.text          = element_text(size = 9),
      legend.title         = element_blank()
    )
  } else {
    theme(legend.position = "none")
  }
  
  p_out <- p +
    labs(title = NULL) +
    base_theme +
    legend_theme
  
  if (!is.null(file)) {
    ggsave(file, p_out,
           width  = width,
           height = height,
           units  = "in",
           dpi    = dpi,
           bg     = "white")
  }
  
  invisible(p_out)
}

#── example ─────────────────────────────────────────────────────────
# rand_plot_clean <- style_eval_plot(plot_rand_dis,
#                                    show_legend = TRUE,
#                                    file = "rand_vs_duration.png")


#── example usage ───────────────────────────────────────────────────────────
rand_plot_clean <- style_eval_plot(plot_rand_dis,
                                   file = "rand_vs_duration_100_real_seed_123.png")

ari_plot_clean  <- style_eval_plot(plot_adj_rand_dis,
                                   file = "ari_vs_duration_100_real_seed_123.png")

#Try avg metric
full_sim_avg <- method_sim %>% evaluate(list(rand_index_metric, adjusted_rand_index_metric))
rand_plot_dis_n100 <- plot_eval_by(full_sim_avg, "rand_index", varying = "disappearance_duration", main = "Cluster Disappearance Simulation", 
             xlab = "Disappearance Duration", ylab = "Rand Index")
plot_eval_by(full_sim_avg, "adjusted_rand_index", varying = "disappearance_duration", main = "Cluster Disappearance Simulation", 
             xlab = "Disappearance Duration", ylab = "Rand Index")
rand_plot_dis_n100_clean <- style_eval_plot(rand_plot_dis_n100,
                                   file = "rand_vs_duration_n100.png")

###############  intersection simulation ###################
hSigma <- 20
n_runs <- 100
ncores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
intersect_clusters <- seq(0.9, 2.5, by = 0.1) #change to 0.5 to 2.5

sim_intersect <- new_simulation(name = "intersection_study",
                      label = "Intersection Study Simulation") %>%
  generate_model(make_intersection_model,
                 intersect_cluster = as.list(intersect_clusters),
                 vary_along = "intersect_cluster") %>%
  simulate_from_model(nsim = n_runs)


method_sim_intersect <- sim_intersect %>% run_method(
  list(kernel_em_CV_method, const_method, hungarian_method),
  parallel = list(
    socket_names = rep("localhost", ncores),    # Number of cores to use
    libraries = c("flowkernel", "foreach", "doParallel", "reshape2", "mclust")   # Libraries to load on each worker
  )
)

full_sim_intersect <- method_sim_intersect %>% evaluate(list(combined_concatenated_rand_index_metric, combined_concatenated_ARI_metric))

## 1) one all-in-one workspace file (easy to reload)
save(sim_intersect, method_sim_intersect, full_sim_intersect,
     file = "rdata/intersect_workspace.RData",
     compress = "xz")          # smallest size, still fast

## 2) individual RDS files (handy if only one object is needed later)
saveRDS(sim_intersect,        file = "rdata/sim_intersect.rds",        compress = "xz")
saveRDS(method_sim_intersect, file = "rdata/method_sim_intersect.rds", compress = "xz")
saveRDS(full_sim_intersect,   file = "rdata/full_sim_intersect.rds",   compress = "xz")


plot_rand_inter <- plot_eval_by(full_sim_intersect, "combined_concatenated_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation", 
             xlab = "Cluster Intersection", ylab = "Rand Index")
plot_adj_rand_inter <- plot_eval_by(full_sim_intersect, "combined_concatenated_adjusted_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation", 
                                xlab = "Cluster Intersection", ylab = "Adjusted Rand Index")

#Try avg metric
full_sim_intersect_avg <- method_sim_intersect %>% evaluate(list(rand_index_metric, adjusted_rand_index_metric))

saveRDS(full_sim_intersect_avg,        file = "rdata/full_sim_intersect_avg.rds",        compress = "xz")

rand_plot_inter_n100 <- plot_eval_by(full_sim_intersect_avg, "rand_index", varying = "intersect_cluster", main = "Cluster Disappearance Simulation", 
                                   xlab = "Cluster Intersection", ylab = "Rand Index")
plot_eval_by(full_sim_intersect_avg, "adjusted_rand_index", varying = "intersect_cluster", main = "Cluster Disappearance Simulation", 
             xlab = "Cluster Intersection", ylab = "Adjusted Rand Index")

library(ggplot2)

## remove every layer whose geom is the horizontal‐error-bar geom
rand_plot_inter_n100$layers[[3]]$aes_params$width  <- 0

library(ggplot2)

## 1. change x-axis limits to 1.5 … 2.5
rand_plot_inter_n100 <- rand_plot_inter_n100 +
  scale_x_continuous(limits = c(1.5, 2.5))

## 2. drop horizontal whiskers (keep vertical bars only)
rand_plot_inter_n100$layers[[3]]$aes_params$width <- 0

## 3. clean & save with your styling helper
rand_plot_inter_n100_clean <- style_eval_plot(
  rand_plot_inter_n100,
  file = "rand_vs_inter_n100_cut.png"
)


rand_plot_inter_n100                # now only vertical error bars remain

rand_plot_inter_n100_clean <- style_eval_plot(rand_plot_inter_n100,
                                            file = "rand_vs_inter_n100_cut.png")

library(ggplot2)

make_vertical_errorbars <- function(p) {
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomErrorbar")) {
      p$layers[[i]]$geom <- GeomLinerange        # swap geom object
      p$layers[[i]]$geom_params$width <- NULL    # no need for width
    }
  }
  p
}

# ── use it ──────────────────────────────────────────────────────────
plot_adj_rand_inter_v <- make_vertical_errorbars(plot_adj_rand_inter)
plot_adj_rand_inter_v




library(ggplot2)

make_vertical_errorbars <- function(p,
                                    show_legend = TRUE,
                                    file        = NULL,
                                    width  = 4,
                                    height = 2.4,
                                    dpi    = 300) {
  
  ## ── 1. drop horizontal caps on every GeomErrorbar layer ───────────
  p2 <- p
  p2$layers <- lapply(p2$layers, function(lyr) {
    if (inherits(lyr$geom, "GeomErrorbar"))
      lyr$geom_params$width <- 0   # <- no caps
    lyr
  })
  
  ## ── 2. reuse the styling helper from before ───────────────────────
  p_final <- style_eval_plot(p2,
                             show_legend = show_legend,
                             file        = file,
                             width       = width,
                             height      = height,
                             dpi         = dpi)
  
  invisible(p_final)
}

# ── examples ──────────────────────────────────────────────────────────
rand_v  <- make_vertical_errorbars(plot_rand_inter,
                                   file = "rand_intersection.png")

ari_v   <- make_vertical_errorbars(plot_adj_rand_inter,
                                   file = "ari_intersection.png",
                                   show_legend = FALSE)   # hide legend if desired

plot_rand_inter_clean <- plot_rand_inter +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
plot_rand_inter_clean

ggsave("Figures/rand_index_comparison_inter.png", plot_rand_inter_clean,
       width = 10, height = 4, dpi = 300)

# strip the caps from every geom_errorbar layer in the plot
plot_rand_inter_clean$layers <- lapply(plot_rand_inter_clean$layers, function(layer) {
  if (inherits(layer$geom, "GeomErrorbar")) layer$aes_params$width <- 0
  layer
})

library(cowplot)

## take one legend
leg <- get_legend(plot_rand_dis_clean)   # legend is identical for both plots

## drop the legends from the two panels
p1 <- plot_rand_dis_clean  + theme(legend.position = "none")
p2 <- plot_rand_inter_clean + theme(legend.position = "none")

## put panels side-by-side and add panel labels
panels  <- plot_grid(
  p1, p2,
  labels = c("(a) Disappearance simulation", "(b) Intersect simulation"),
  label_size = 14,
  nrow = 1, align = "h"
)

## attach the single legend to the right
final_fig <- plot_grid(panels, leg, ncol = 2, rel_widths = c(1, .22))

ggsave("Figures/rand_index_comparison_dis.png", final_fig,
       width = 10, height = 4, dpi = 300)


full_sim_intersect <- method_sim_intersect %>% evaluate(list(rand_index_metric, adjusted_rand_index_metric, combined_concatenated_rand_index_metric))


full_sim_intersect <- method_sim_intersect %>% evaluate(list(combined_aligned_label_error_metric))

# Save both simulation objects together
save(sim_intersect, method_sim_intersect, full_sim_intersect, file = "inter_n3_objects.RData")

plot_eval_by(full_sim_intersect, "adjusted_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "combined_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")

plot_eval_by(full_sim_intersect, "combined_aligned_label_error", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")


p <- plot_eval_by(
  full_sim_intersect, 
  metric_name = "combined_concatenated_rand_index", 
  varying     = "intersect_cluster",
  type        = "aggregated",
  main        = "Intersecting Cluster Simulation", 
  xlab        = "Cluster Intersection", 
  ylab        = "Rand Index",
  use_ggplot2 = TRUE
)

# Layer 3 is the error bars. Set their width to 0 to remove horizontal whiskers:
p$layers[[3]]$geom_params$width <- 0
p$layers[[3]]$aes_params$width  <- 0


# Now reprint or save p
print(p)


#Debug

met <- kernel_em_CV_method  # "met" is now your method object
met <- hungarian_method

m <- model(sim_intersect, subset = "intersection_model_intersect_1.8/intersect_cluster_1.8")
d <- draws(sim_intersect, subset = "intersection_model_intersect_1.8/intersect_cluster_1.8", index = 1)
.Random.seed <<- as.integer(c(10407, 1926221190, -1343718741, 518015349, -127955897, 1365626998, 1718833432))

# Now run the method code on that one model/draw
met@method(model = m, draw = d@draws$r1.3)



#Try from loaded sim

full_sim_intersect <- sim %>% evaluate(list(rand_index_metric, adjusted_rand_index_metric, combined_rand_index_metric))

full_sim_intersect <- sim %>% evaluate(list(combined_concatenated_rand_index_metric))
full_sim_intersect <- sim %>% evaluate(list(combined_aligned_label_error_metric))

plot_eval_by(sim, "adjusted_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(sim, "rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "combined_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(sim, "combined_concatenated_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
View(sim)

method_sim_intersect <- sim %>% run_method(
  list(kernel_em_CV_method, const_method, hungarian_method),
  parallel = list(
    socket_names = ncores,    # Number of cores to use
    libraries = c("flowkernel", "foreach", "doParallel", "reshape2", "mclust")   # Libraries to load on each worker
  )
)

## Normalized eval plot ## 
library(dplyr)
library(tidyr)
library(ggplot2)

# Convert the evaluation results to a data frame.
eval_df <- as.data.frame(evals(full_sim_intersect))

# Normalize the combined concatenated Rand index within each Model/Draw group.
normalized_df <- eval_df %>%
  group_by(Model, Draw) %>% 
  mutate(norm_combined_rand_index = combined_concatenated_rand_index /
           combined_concatenated_rand_index[Method == "kernel_em_cv"]) %>%
  ungroup()

# Separate the 'Model' column to extract the intersection level.
# This assumes the 'Model' column is in the form
# "intersection_model_intersect_X/intersect_cluster_X"
normalized_df <- normalized_df %>%
  separate(Model, into = c("model", "intersect_cluster"), sep = "/") %>%
  # Remove the prefix and convert to numeric.
  mutate(intersect_cluster = as.numeric(sub("intersect_cluster_", "", intersect_cluster)))

# Aggregate by Method and intersection level: compute the mean and standard error.
agg_df <- normalized_df %>%
  group_by(Method, intersect_cluster) %>%
  summarise(avg_norm_rand = mean(norm_combined_rand_index, na.rm = TRUE),
            se_norm_rand = sd(norm_combined_rand_index, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

# Plot the average normalized metric with error bars.
ggplot(agg_df, aes(x = intersect_cluster, y = avg_norm_rand, color = Method, group = Method)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = avg_norm_rand - se_norm_rand, ymax = avg_norm_rand + se_norm_rand),
                width = 0.05) +
  labs(title = "Normalized Rand Index (relative to kernel_em_cv)",
       x = "Intersection Level",
       y = "Average Normalized Rand Index") +
  scale_x_continuous(breaks = unique(agg_df$intersect_cluster)) +
  theme_minimal()

#Function

normalized_plot <- function(sim_obj, metric, varying = "intersect_cluster", main = "Intersecting Cluster Simulation") {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  
  # Extract the evaluation results as a list and convert to a data frame.
  eval_list <- evals(sim_obj)
  eval_df <- map_df(eval_list, ~ as.data.frame(.x))
  
  # Normalize the specified metric within each simulation draw.
  # We assume that the kernel EM method is labeled "kernel_em_cv".
  normalized_df <- eval_df %>%
    group_by(Model, Draw) %>%
    mutate(norm_metric = .data[[metric]] / .data[[metric]][Method == "kernel_em_cv"]) %>%
    ungroup()
  
  # Extract the intersection level from the 'Model' column.
  # Assumes Model is in the format "intersection_model_intersect_X/intersect_cluster_X"
  normalized_df <- normalized_df %>%
    separate(Model, into = c("model", "intersect_cluster"), sep = "/") %>%
    mutate(intersect_cluster = as.numeric(sub("intersect_cluster_", "", intersect_cluster)))
  
  # Aggregate: compute mean and standard error for each Method at each intersection level.
  agg_df <- normalized_df %>%
    group_by(Method, intersect_cluster) %>%
    summarise(avg_norm = mean(norm_metric, na.rm = TRUE),
              se_norm = sd(norm_metric, na.rm = TRUE) / sqrt(n()),
              .groups = "drop")
  
  # Create the plot.
  p <- ggplot(agg_df, aes(x = intersect_cluster, y = avg_norm, color = Method, group = Method)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = avg_norm - se_norm, ymax = avg_norm + se_norm),
                  width = 0.05) +
    labs(title = main,
         x = "Intersection Level",
         y = paste("Normalized", metric)) +
    scale_x_continuous(breaks = sort(unique(agg_df$intersect_cluster))) +
    theme_minimal()
  
  return(p)
}
# For adjusted Rand index
p1 <- normalized_plot(full_sim_intersect, "adjusted_rand_index", main = "Adjusted Rand Index")
print(p1)

# For Rand index
p2 <- normalized_plot(full_sim_intersect, "rand_index", main = "Rand Index")
print(p2)

# For combined Rand index
p3 <- normalized_plot(full_sim_intersect, "combined_rand_index", main = "Combined Rand Index")
print(p3)


p4 <- normalized_plot(full_sim_intersect, "combined_concatenated_rand_index", main = "Intersecting Cluster Simulation")
print(p4)
plot_eval_by(full_sim_intersect, "combined_concatenated_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")

## ## ## ## ## ## ## ## ## ## ## ## 

cv_results_mid <- list()  # Initialize list to store CV results for each duration

for (dur in disappearance_durations) {
  # Create a temporary simulation to get a single draw
  sim_temp <- new_simulation(name = sprintf("temp_sim_%s", dur),
                             label = sprintf("Temporary Simulation for Dur %s", dur)) %>%
    generate_model(make_disappearance_model,
                   transition_duration = 30,
                   disappearance_duration = dur) %>%
    simulate_from_model(nsim = 1)
  
  # Extract one draw to use for cross-validation
  d <- draws(sim_temp, index = 1)  # Get the draws component
  dat <- d@draws$r1.1  # Access the specific draw in the @draws slot
  
  # Perform cross-validation
  cv_result <- cv_log_grid_search_parallel(data = dat$y, K = 2, hSigma = 15,
                                           dates = dat$dates, biomass = dat$biomass,
                                           leave_out_every = 5, grid_size = 10)
  cv_results_mid[[as.character(dur)]] <- list(hmu = cv_result$best_hmu, hpi = cv_result$best_hpi)
}

cv_results_mid <- cv_results_mid

# New intersection model from o3

#Create intersection simulation object
hSigma <- 20
n_runs <- 3
intersect_clusters <- seq(0.5, 2.5, by = 0.3)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_ON_NODE"))
sim_intersect <- new_simulation(name = "intersection_study",
                                label = "Intersection Study Simulation") %>%
  generate_model(make_intersection_model_flip,
                 intersect_cluster = as.list(intersect_clusters),
                 vary_along = "intersect_cluster") %>%
  simulate_from_model(nsim = n_runs)


method_sim_intersect <- sim_intersect %>% run_method(
  list(kernel_em_CV_method, const_method, hungarian_method),
  parallel = list(
    socket_names = rep("localhost", ncores),    # Number of cores to use
    libraries = c("flowkernel", "foreach", "doParallel", "reshape2", "mclust")   # Libraries to load on each worker
  )
)

full_sim_intersect <- method_sim_intersect %>% evaluate(list(rand_index_metric, adjusted_rand_index_metric, combined_rand_index_metric, combined_concatenated_rand_index_metric))

full_sim_intersect <- method_sim_intersect %>% evaluate(list(combined_concatenated_rand_index_metric))
full_sim_intersect <- method_sim_intersect %>% evaluate(list(combined_aligned_label_error_metric))

# Save both simulation objects together
save(sim_intersect, method_sim_intersect, full_sim_intersect, file = "inter_flip_n10_objects.RData")

plot_eval_by(full_sim_intersect, "adjusted_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "combined_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "combined_concatenated_rand_index", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")
plot_eval_by(full_sim_intersect, "combined_aligned_label_error", varying = "intersect_cluster", main = "Intersecting Cluster Simulation")


