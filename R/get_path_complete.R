#' @title Get Complete Paths from a Dendrogram
#'
#' @description
#' This function computes recommended paths from a dendrogram structure using parallel processing. 
#' It utilizes the `find_recommended_path` function to identify optimal paths based on a minimum 
#' cardinality constraint, distributing the computation across multiple worker nodes.
#'
#' @param minidend A dendrogram structure from which to derive paths. This should be created from 
#'                 a hierarchical clustering result.
#' @param window_data A data frame or matrix containing the data associated with the nodes in the dendrogram. 
#'                    This data is used for path recommendations.
#' @param min_card An integer specifying the minimum number of leaves (or nodes) that must be present 
#'                  in a path for it to be considered valid.
#' @param worker_number An integer representing the number of worker nodes to be used for parallel processing.
#'
#' @return A list where each element contains the recommended paths for the corresponding node in the dendrogram. 
#'         Each path includes information about the nodes and their associated scores.
#'
#' @details
#' The function creates a cluster of worker nodes, loads necessary libraries, and exports required 
#' variables and functions to each worker. It then applies the `find_recommended_path` function 
#' in parallel to the leaves of the provided dendrogram, gathering results into a single list. 
#' Finally, the cluster is stopped, and the results are returned.
#'
#' @export
get_path_complete <- function(minidend, window_data, min_card,worker_number){

  # Set up a cluster
  cl <- makeCluster(worker_number)
  
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