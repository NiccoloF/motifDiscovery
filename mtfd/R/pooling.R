#' @title probKMA_plot
#'
#' @description Plot the results of probKMA.
#' @param dataMatrix matrix that must be extended  
#' @param maxLen extension length
#' @return matrix extended up to maxLen
#' @author Niccolò Feresini 
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