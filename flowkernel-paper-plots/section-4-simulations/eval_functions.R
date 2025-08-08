## @knitr metrics

# Define Rand Index metric
rand_index_metric <- new_metric(name = "rand_index",
                                label = "Rand Index",
                                metric = function(model, out) {
                                  true_z <- out$true_z
                                  fit_resp <- out$fit$resp
                                  calculate_mean_rand_index(true_z, fit_resp)
                                })

# Define Adjusted Rand Index metric
adjusted_rand_index_metric <- new_metric(name = "adjusted_rand_index",
                                         label = "Adjusted Rand Index",
                                         metric = function(model, out) {
                                           true_z <- out$true_z
                                           fit_resp <- out$fit$resp
                                           calculate_adj_mean_rand_index(true_z, fit_resp)
                                         })

combined_rand_index_metric <- new_metric(name = "combined_rand_index",
                                label = "Combined Rand Index",
                                metric = function(model, out) {
                                  true_z <- out$true_z
                                  fit_resp <- out$fit$resp
                                  combined_rand_index(true_z, fit_resp)
                                })



library(fossil)
library(mclust)

sample_resp <- function(prob_row) {
  n <- length(prob_row)
  sampled_index <- sample(1:n, size = 1, prob = prob_row)
  one_hot <- numeric(n)
  one_hot[sampled_index] <- 1
  return(one_hot)
}

# Function to convert a response matrix to a one-hot encoded matrix
z_mat <- function(resp_matrix) {
  result_matrix <- t(apply(resp_matrix, 1, sample_resp))
  return(result_matrix)
}

# Function to convert one-hot encoded matrices to cluster assignments
convert_to_assignments <- function(one_hot_matrix) {
  assignments <- apply(one_hot_matrix, 1, which.max)
  return(assignments)
}

# Function to compute the mean Rand index for lists of true and predicted assignments
calculate_mean_rand_index <- function(true_list, resp_list) {
  # Convert response matrices to one-hot and then to assignments
  one_hot_list <- lapply(resp_list, z_mat)
  assignments_list <- lapply(one_hot_list, convert_to_assignments)
  
  # Calculate the Rand index for each pair and compute the mean
  rand_ind_sum <- 0
  n <- length(true_list)
  for (i in 1:n) {
    rand_ind <- rand.index(true_list[[i]], assignments_list[[i]])
    #print(rand_ind)  # Optional: print each Rand index
    rand_ind_sum <- rand_ind_sum + rand_ind
  }
  
  rand_ind_mean <- rand_ind_sum / n
  return(rand_ind_mean)
}


#Adjusted rand index

calculate_adj_mean_rand_index <- function(true_list, resp_list) {
  # Convert response matrices to one-hot and then to assignments
  one_hot_list <- lapply(resp_list, z_mat)
  assignments_list <- lapply(one_hot_list, convert_to_assignments)
  
  # Calculate the Rand index for each pair and compute the mean
  adj_rand_ind_sum <- 0
  n <- length(true_list)
  for (i in 1:n) {
    adj_rand_ind <- adjustedRandIndex(true_list[[i]], assignments_list[[i]])
    #print(adj_rand_ind)  # Optional: print each Rand index
    adj_rand_ind_sum <- adj_rand_ind_sum + adj_rand_ind
  }
  
  adj_rand_ind_mean <- adj_rand_ind_sum / n
  return(adj_rand_ind_mean)
}

# Function to calculate the combined Rand index
combined_rand_index <- function(true_list, resp_list) {
  # Convert response matrices to one-hot encoded matrices
  one_hot_list <- lapply(resp_list, z_mat)
  
  # Convert one-hot encoded matrices to cluster assignments
  assignments_list <- lapply(one_hot_list, convert_to_assignments)
  
  # Concatenate all true assignments and predicted assignments
  true_assignments <- unlist(true_list)
  predicted_assignments <- unlist(assignments_list)
  
  # Compute the Rand index between the concatenated assignments
  rand_ind <- rand.index(true_assignments, predicted_assignments)
  
  return(rand_ind)
}


#New metrix Feb 2025

#Concatenate and use Rand
combined_concatenated_rand_index_metric <- new_metric(
  name = "combined_concatenated_rand_index",
  label = "Combined Concatenated Rand Index",
  metric = function(model, out) {
    # Concatenate the true labels from all time points.
    true_assignments <- unlist(out$true_z)
    
    # For each time point, get the predicted assignments using max.col,
    # then concatenate the results.
    predicted_assignments <- unlist(lapply(out$fit$resp, function(resp_matrix) {
      max.col(resp_matrix)
    }))
    
    # Compute the Rand index between the concatenated true and predicted labels.
    rand_index_value <- rand.index(true_assignments, predicted_assignments)
    
    return(rand_index_value)
  }
)

# Combined, concatenated **Adjusted** Rand Index
combined_concatenated_ARI_metric <- new_metric(
  name  = "combined_concatenated_adjusted_rand_index",
  label = "Combined Concatenated Adjusted Rand Index",
  metric = function(model, out) {
    ## concatenate true labels across all time points
    true_assignments <- unlist(out$true_z)
    
    ## concatenate predicted labels (MAP, via max.col on each resp matrix)
    predicted_assignments <- unlist(
      lapply(out$fit$resp, function(resp_matrix) max.col(resp_matrix))
    )
    
    ## adjusted Rand index
    adjustedRandIndex(true_assignments, predicted_assignments)
  }
)


# Unsued metric:
# 
# # install.packages("clue")
# library(clue)
# 
# # This function aligns predicted labels to true labels and returns the misclassification error.
# global_aligned_error <- function(true_labels, pred_labels) {
#   # Ensure that labels are factors
#   true_labels <- as.factor(true_labels)
#   pred_labels <- as.factor(pred_labels)
#   
#   # Build the contingency table between predicted and true labels.
#   cont_matrix <- table(pred_labels, true_labels)
#   
#   # For the assignment problem, we want a square cost matrix.
#   K <- max(nrow(cont_matrix), ncol(cont_matrix))
#   
#   # Pad the contingency matrix with zeros if necessary.
#   if (nrow(cont_matrix) < K) {
#     cont_matrix <- rbind(cont_matrix, matrix(0, nrow = K - nrow(cont_matrix), ncol = ncol(cont_matrix)))
#   }
#   if (ncol(cont_matrix) < K) {
#     cont_matrix <- cbind(cont_matrix, matrix(0, nrow = nrow(cont_matrix), ncol = K - ncol(cont_matrix)))
#   }
#   
#   # Convert the maximization problem into a minimization problem:
#   cost_matrix <- max(cont_matrix) - cont_matrix
#   
#   # Solve the assignment problem.
#   assignment <- solve_LSAP(cost_matrix)
#   
#   # Build a mapping from predicted labels to true labels.
#   pred_levels <- levels(pred_labels)
#   mapping <- setNames(as.integer(assignment[seq_along(pred_levels)]), pred_levels)
#   
#   # Reassign the predicted labels using the mapping.
#   aligned_pred_labels <- as.integer(mapping[as.character(pred_labels)])
#   
#   # Convert true labels to integers.
#   true_labels_int <- as.integer(true_labels)
#   
#   # Compute the misclassification error.
#   error_rate <- mean(aligned_pred_labels != true_labels_int)
#   
#   return(error_rate)
# }
# 
# # Now define a new metric that uses this function over all time points.
# # This metric first concatenates the assignments (from max.col of each responsibility matrix)
# # and then computes the global aligned error.
# combined_aligned_label_error_metric <- new_metric(
#   name = "combined_aligned_label_error",
#   label = "Combined Aligned Label Error",
#   metric = function(model, out) {
#     # Get the true assignments across time (assumed to be provided as a list, one per time point)
#     true_assignments <- unlist(out$true_z)
#     
#     # Get the predicted assignments.
#     # Here we assume that for each time point, the responsibilities matrix has one row per observation.
#     predicted_assignments <- unlist(lapply(out$fit$resp, function(resp_matrix) {
#       # Use max.col to choose the cluster with the highest probability.
#       max.col(resp_matrix)
#     }))
#     
#     # Compute the globally aligned error.
#     error_rate <- global_aligned_error(true_assignments, predicted_assignments)
#     
#     return(1 - error_rate)
#   }
# )





