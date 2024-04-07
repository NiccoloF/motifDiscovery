#' @title get_accolites
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

get_accolites <- function(leaf_label, window_data, portion_len, multiple){
  # number of overlapping elements that define accolites
  overlap <- floor(portion_len/2)
  leaf_label_curve <-  (strsplit(leaf_label, '_') %>%
                          unlist() %>%
                          as.numeric())[1]
  
  leaf_index  <- which(rownames(window_data) == leaf_label) # index in window_data
  r_accolites <- (leaf_index+1):min(leaf_index+overlap, dim(window_data)[1])  # accolites to the right
  l_accolites <- (leaf_index-1):max(leaf_index-overlap,0)  # accolites to the left
  accolites   <- c(leaf_index, l_accolites, r_accolites) # overall accolites
  
  # get the accolites
  leaf_accolites <- rownames(window_data)[accolites]
  leaf_accolites <- leaf_accolites[!is.na(leaf_accolites)] #remove NAs
  # check if they come from the same curve (multiple curves cases)
  if(multiple == TRUE){
    leaf_accolites_curve <- strsplit(leaf_accolites, '_') %>%
      unlist() %>%
      as.numeric() %>%
      matrix(ncol=3, byrow=T) %>%
      as.data.frame()
    to_delete <- which(leaf_accolites_curve$V1 != leaf_label_curve)
    if(length(to_delete) > 0){
      leaf_accolites <- leaf_accolites[-to_delete]
    }
  }
  return(leaf_accolites)
}