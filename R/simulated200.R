#' @title Simulated data for local clustering 
#'
#' @description
#' 18 curves of length 200 (corresponding to 201 evaluation points), each
#' containing exactly one motif of length 60 (corresponding to 61 evaluation
#' points), with noise level sigma=0.1
#' Curves 1-6, 14-16 contain one occurrence of motif 1
#' Curves 7-13, 17-18 contain one occurrence of motif 2
#' The data have been generated using a B-spline basis of order 3 and knots at
#' distance 10, following the simulation procedure presented in Cremona and
#' Chiaromonte (Comparison with non-sparse and sparse functional clustering
#' simulation section).
#' @name simulated200
#' @usage data(simulated200)
"simulated200"