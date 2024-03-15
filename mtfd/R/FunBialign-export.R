#' @docType class
#' @title FunBialign class
#' 
#' @author Niccolo' Feresini 
#' @name FunBialign
#' @export FunBialign

# Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
loadModule(module = "FunBialignModule", TRUE)