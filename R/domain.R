#' @title Domain Length Calculation for Curves
#' 
#' @description
#' Computes the length of the domain of a curve by checking the presence of non-NA values in the specified matrix.
#' The function returns a logical vector indicating which rows contain at least one non-NA element, allowing for an evaluation of the effective length of the domain.
#' 
#' @param v A list containing two elements:
#'   \itemize{
#'     \item \code{v[[1]]}: A matrix representing the values of the curve (v(x)).
#'     \item \code{v[[2]]}: A matrix representing the values of the derivative of the curve (v'(x)).
#'   }
#'   Each matrix should have \code{d} columns, where \code{d} is the dimensionality of the curve.
#' 
#' @param use0 A logical value indicating which element of \code{v} to consider:
#'   \itemize{
#'     \item If \code{use0 = TRUE}, the function considers \code{v[[1]]} (the curve values).
#'     \item If \code{use0 = FALSE}, the function considers \code{v[[2]]} (the derivative values).
#'   }
#' 
#' @return A logical vector where each element is \code{TRUE} if at least one element in the corresponding row of the selected matrix (either \code{v[[1]]} or \code{v[[2]]}) is not NA, and \code{FALSE} otherwise.
#' 
#' @details
#' The function evaluates each row of the specified matrix (either \code{v[[1]]} or \code{v[[2]]}) to determine if there is at least one non-NA entry. This is particularly useful for assessing the effective domain of curves where some values may be missing.
#' 
#' @export
.domain <- function(v,use0){
  if(use0){
    rowSums(!is.na(v[[1]]))!=0
  }else{
    rowSums(!is.na(v[[2]]))!=0
  }
}