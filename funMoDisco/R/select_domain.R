#' Select Domain of a Motif
#'
#' @title .select_domain
#'
#' @description
#' This function selects the portion of a motif that is free from NA values based on a specified domain. 
#' It allows for the selection of two matrices, `v0` and `v1`, depending on the boolean flags provided.
#'
#' @param v A list containing two elements: 
#'          \itemize{
#'            \item v[[1]]: A matrix representing \( v(x) \) with \( d \) columns.
#'            \item v[[2]]: A matrix representing \( v'(x) \) with \( d \) columns.
#'          }
#' @param v_dom A boolean vector indicating the domain for `v0` and `v1`. 
#' @param use0 A boolean value indicating whether to select `v0`. If TRUE, `v0` is selected based on `v_dom`.
#' @param use1 A boolean value indicating whether to select `v1`. If TRUE, `v1` is selected based on `v_dom`.
#'
#' @return A list containing the selected portions of `v0` and `v1` (if applicable), with dimensions adjusted accordingly.
#'
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#'
#' @export
.select_domain <- function(v,v_dom,use0,use1){
  if(use0)
    v[[1]]=as.matrix(v[[1]][v_dom,])
  if(use1)
    v[[2]]=as.matrix(v[[2]][v_dom,])
  return(v)
}