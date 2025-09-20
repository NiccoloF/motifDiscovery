#' @title Pad a Matrix to a Specified Number of Rows
#'
#' @description 
#' This function pads a matrix with `NA` values to ensure that the matrix reaches a specified number of rows (`maxLen`). 
#' If the matrix has fewer rows than the desired length, additional rows of `NA` are appended to the matrix. 
#' Otherwise, the matrix is returned unchanged.
#'
#' @param dataMatrix A matrix to be padded. Each row represents a data point.
#' @param maxLen An integer specifying the desired number of rows in the matrix.
#'
#' @details 
#' The function compares the number of rows in `dataMatrix` with `maxLen`. If the matrix has fewer rows, 
#' it pads the matrix with `NA` values until it reaches the specified number of rows. 
#' No changes are made if `dataMatrix` has rows equal to or greater than `maxLen`.
#'
#' @return A matrix padded with `NA` values up to the specified number of rows. 
#' If `dataMatrix` already has the desired number of rows or more, the original matrix is returned.
#'
#' @export
padding <- function(dataMatrix, maxLen) 
{
  actual_len <- nrow(dataMatrix)
  if (actual_len < maxLen) 
  {
    extension_len <- maxLen - actual_len
    extension <- matrix(NA, nrow = extension_len, ncol = ncol(dataMatrix))
    full_mat <- rbind(dataMatrix, extension)
    return(full_mat)
  }else{return(dataMatrix)}
}
