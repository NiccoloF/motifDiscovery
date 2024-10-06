#' @title Custom mapply Function for Parallel Processing
#'
#' @description
#' A wrapper function that conditionally applies either the classical `mapply` function 
#' or `parallel::clusterMap` based on whether a cluster object is provided. 
#' If the number of clusters is not null, it applies `clusterMap` for parallel execution; 
#' otherwise, it defaults to the classical `mapply`.
#'
#' @param cl A cluster object created by `makeCluster`. If `NULL`, the classical `mapply` function is used.
#' @param FUN A function to apply to the arguments.
#' @param ... Arguments to be passed to `FUN`. These should be vectors of equal length.
#' @param MoreArgs A list of additional arguments to be passed to `FUN`.
#' @param SIMPLIFY A logical value indicating whether to simplify the result if possible. Default is `TRUE`.
#' @param USE.NAMES A logical value indicating whether to use names from the first argument. Default is `TRUE`.
#'
#' @return A vector or list containing the results of applying the function `FUN` to the provided arguments. 
#'         The output type depends on the value of `SIMPLIFY`.
#'
#' @details
#' This function is useful for switching between parallel and non-parallel execution 
#' based on the availability of a cluster, allowing for more flexible and efficient code.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' library(parallel)
#' cl <- makeCluster(4) # Create a cluster with 4 workers
#'
#' # Define a simple function to square numbers
#' square_func <- function(x) {
#'   return(x^2)
#' }
#'
#' # Using .mapply_custom with a cluster
#' result_parallel <- .mapply_custom(cl, square_func, 1:5, MoreArgs = list())
#' print(result_parallel) # Prints squared values
#'
#' # Using .mapply_custom without a cluster
#' result_serial <- .mapply_custom(NULL, square_func, 1:5, MoreArgs = list())
#' print(result_serial) # Prints squared values
#'
#' # Stop the cluster after use
#' stopCluster(cl)
#' }
#'
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#'
#' @export
.mapply_custom <- function(cl,FUN,...,MoreArgs=NULL,SIMPLIFY=TRUE,USE.NAMES=TRUE){
  if(is.null(cl)){
    mapply(FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }else{
    parallel::clusterMap(cl,FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }
}