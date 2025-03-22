#' @title Generate Minimum Dendrogram from Hierarchical Clustering
#'
#' @description
#' This function generates a minimum dendrogram by performing hierarchical clustering on an input distance matrix. 
#' The function identifies an optimal height cut based on the largest gap in cluster heights and returns the resulting 
#' dendrogram, allowing for further analysis of clustering structures.
#'
#' @param adj_fMSR A numeric matrix or a distance object representing the dissimilarity matrix 
#'                  computed from functional data. This matrix is used to perform hierarchical clustering.
#'
#' @return A dendrogram object representing the clustered data. The dendrogram is cut at the identified height,
#'         resulting in a tree structure that can be used for further analysis or visualization.
#'
#' @details
#' The function performs the following steps:
#' 1. Uses `fastcluster::hclust` to perform hierarchical clustering on the provided dissimilarity matrix using 
#'    the "complete" method.
#' 2. Extracts the heights of the clusters and identifies the largest gap in heights to determine the optimal 
#'    cut position for the dendrogram.
#' 3. Cuts the dendrogram at the identified height and returns the lower part of the cut.
#'
#' @importFrom fastcluster hclust
#' @export
get_minidend <- function(adj_fMSR){
  # generate dendrogram
  total_hc <- fastcluster::hclust(adj_fMSR, method = "complete")
  total_dend <- as.dendrogram(total_hc) # transform in dendrogram
  
  # Identify height_cut
  # It is made by ordering all the heights and identifying the largest gap.
  # The height cut corresponds to the mean between the elements generating the gap
  cut_position <- sort(total_hc$height) %>% diff() %>% which.max() # total_hc$height is filled with the cluster distance for each step
  height_cut <- sort(total_hc$height)[c(cut_position, cut_position + 1)] %>%
    mean()
  
  # Generating minidend
  minidend <- cut(total_dend, h = height_cut)$lower
  return(minidend)
}