#' @title get_path_complete
#'
#' @description Run multiple times probKMA function with different K,c and initializations,
#' with the aim to find a set of candidate motifs.
#' If the folder name_KK_cc is already present and n result files are already present,
#' load them and continue with the n_init-n runs.
#'
#' @param Y0 list of N vectors, for univariate curves y_i(x), or list of N matrices with d columns,
#' for d-dimensional curves y_i(x),  with the evaluation of curves (all curves should be evaluated
#' on a uniform grid). When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
#' @param Y1 list of N vectors, for univariate derivative curves y'_i(x), or
#' list of N matrices with d columns, for d-dimensional derivatibe curves y'_i(x),
#' with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
#' When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions.
#' Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
#' @param K vector with numbers of motifs that must be tested.
#' @param c vector with minimum motifs lengths that must be tested.
#' @param n_init number of random initialization for each combination of K and c.
#' @param name name of the folders when the results are saved.
#' @param names_var vector of length d, with names of the variables in the different dimensions.
#' @param probKMA_options list with options for probKMA (see the help of probKMA).
#' @param silhouette_align True or False. If True, try all possible alignments between the curve pieces
#' when calculating the adapted silhouette index on the results of probKMA
#' @param plot if TRUE, summary plots are drawn.
#' @return A list containing: K, c, n_init and name;...
#' @return \item{times}{ list of execution times of ProbKMA for each combination of K, c, and n_init}
#' @return \item{silhouette_average_sd}{ list of the mean (silhouette_average) and standard deviation (silhouette_sd) of the silhouette indices for each execution of the ProbKMA function}
#' @author Marzia Angela Cremona & Francesca Chiaromonte

get_path_complete <- function(minidend, window_data, min_card){
  # Determine the number of cores
  num_cores <- detectCores() - 1
  # Set up a cluster
  cl <- makeCluster(num_cores)
  
  # Load necessary libraries on each worker
  clusterEvalQ(cl, {
    library(dplyr) 
    library(data.table)
  })
  
  # Export necessary variables and functions to the cluster
  clusterExport(cl, c("find_recommended_path", "window_data", "min_card","partition_leaves", "get_nodes_xy")
                ,envir = environment())
  
  # Use parLapply to run in parallel
  all_paths <- parLapply(cl, minidend, find_recommended_path, 
                         window_data = window_data, min_card = min_card)
  
  # Stop the cluster
  stopCluster(cl)
  
  return(all_paths)
}