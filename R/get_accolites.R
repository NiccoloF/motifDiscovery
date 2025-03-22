#' @title Get Accolites for a Given Leaf Label
#'
#' @description
#' This function retrieves the accolites (adjacent elements) of a specified leaf label from a dataset. 
#' Accolites are defined as the elements that overlap with the specified portion length around the leaf.
#' The function can filter accolites based on their origin curve in cases where multiple curves are present.
#'
#' @param leaf_label A character string representing the label of the leaf for which to find accolites.
#' @param window_data A data frame or matrix where the row names correspond to the labels of elements,
#'                    containing the data from which accolites will be extracted.
#' @param portion_len An integer specifying the total length of the portion used to define the accolites. 
#'                    The overlap is computed as half of this length.
#' @param multiple A logical indicating whether to check for multiple curves. If TRUE, the function 
#'                 filters out accolites that do not originate from the same curve as the specified leaf label.
#'
#' @return A character vector containing the labels of the identified accolites. 
#'         If no accolites are found, the function returns an empty vector.
#'
#' @details
#' The function works as follows:
#' 1. It calculates the index of the specified leaf label within the provided `window_data`.
#' 2. Determines the range of indices representing the accolites by calculating the overlap based on `portion_len`.
#' 3. Retrieves the corresponding leaf labels from the `window_data`.
#' 4. If the `multiple` argument is TRUE, it checks if the accolites come from the same curve as the leaf label, 
#'    removing those that do not.
#'
#' @export
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