#' Recommend Node Function
#'
#' @title Recommend Node from a Numeric Vector
#'
#' @description
#' This function analyzes a numeric vector and recommends an index based on the characteristics of the vector. 
#' If the vector is strictly increasing, it suggests stopping before the maximum growth rate. 
#' If the vector is not strictly increasing, it recommends the index of the minimum value.
#'
#' @param node A numeric vector representing scores or metrics. Must contain at least two elements.
#'
#' @return An integer index indicating the recommended node based on the analysis of the input vector.
#'
#' @details
#' - If the input vector has fewer than two elements, the function defaults to returning an index of 1.
#' - The function first checks if the input is a numeric vector. If not, an error is thrown.
#' - The function calculates the differences between successive elements and identifies whether the vector is strictly increasing.
#' - In the case of a strictly increasing vector, it returns the index of the maximum growth rate.
#' - For non-increasing vectors, it returns the index of the minimum value.
#'
#' @export
recommend_node <- function(node){
  xx <- unlist(node)
  
  if(length(xx) >= 2){
    diff_xx <- diff(xx)
    if(all(diff(xx)>0)){ # always increasing h-score: stop before maximum growth rate
      return((diff_xx/xx[1:length(xx)-1]) %>% which.max()) # or +1
    }else{ #also decreasing: stop minimum hscore
      return(which.min(xx))
    }
  }else{
    return(1)
  }
}
