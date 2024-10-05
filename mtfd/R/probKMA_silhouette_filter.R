#' @title Filter Motifs from probKMA Results Based on Silhouette and Size Thresholds
#'
#' @description This function filters the motifs identified by the `probKMA` algorithm based on a threshold for the average silhouette index and a minimum size criterion. Motifs that meet or exceed both the silhouette index threshold and the minimum number of curves (size) are retained. The function returns a cleaned version of the input data, including the filtered motifs, their derivatives, and associated matrices.
#'
#' @param probKMA_results A list representing the output of the `probKMA` function with `return_options = TRUE`. This output includes the motifs, dissimilarity matrices, and membership matrices required for filtering.
#' @param silhouette A list output from the `probKMA_silhouette` function, which provides silhouette scores for each motif.
#' @param sil_threshold A numeric value representing the threshold for the average silhouette index. Motifs with a silhouette index greater than or equal to this value will be retained. Default is 0.5.
#' @param size_threshold An integer representing the minimum number of curves (size) in a motif cluster. Motifs with a cluster size greater than or equal to this value will be retained. Default is 2.
#' 
#' @return A list containing the filtered results:
#' \item{V0_clean}{Filtered motifs.}
#' \item{V1_clean}{Derivatives of the filtered motifs.}
#' \item{D}{Dissimilarity matrix.}
#' \item{D_clean}{Filtered dissimilarity matrix after cleaning motifs.}
#' \item{P}{Membership matrix.}
#' \item{P_clean}{Filtered membership matrix after cleaning motifs.}
#' \item{c}{Filtered minimum motif lengths.}
#' \item{K}{Vector containing the number of motifs repeated by the filtered motifs.}
#'
#' @export
probKMA_silhouette_filter <- function(probKMA_results,silhouette,sil_threshold=0.5,size_threshold=2){
  index_sil=which(silhouette$silhouette_average>=sil_threshold)
  index_size=which(colSums(probKMA_results$P_clean)>=size_threshold)
  index=intersect(index_sil,index_size)
  if(length(index)==0)
    return(NULL)
  return(list(V0_clean=probKMA_results$V0_clean[index],V1_clean=probKMA_results$V1_clean[index],
              D=probKMA_results$D[,index],D_clean=probKMA_results$D_clean[,index],P=probKMA_results$P[,index],P_clean=probKMA_results$P_clean[,index],
              c=probKMA_results$c[index],K=rep(probKMA_results$K,length(index))))
}