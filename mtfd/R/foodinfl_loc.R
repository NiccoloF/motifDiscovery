#' Simulated database for local clustering 
#'
#' 18 curves of length 200 (corresponding to 201 evaluation points), each
#' containing exactly one motif of length 60 (corresponding to 61 evaluation
#' points), with noise level sigma=0.1
#' Curves 1-6, 14-16 contain one occurrence of motif 1
#' Curves 7-13, 17-18 contain one occurrence of motif 2
#' The data have been generated using a B-spline basis of order 3 and knots at
#' distance 10, following the simulation procedure presented in Cremona and
#' Chiaromonte (Comparison with non-sparse and sparse functional clustering
#' simulation section).
#'
#' @format List with the 10 following elements:
#' \describe{
#' \item{\code{Y0}}{list of 18 vectors with the evaluation of the 18 curves}
#' \item{\code{Y1}}{list of 18 vectors with the evaluation of the 18 curve
#' derivatives}
#' \item{\code{dist_knots}}{distance between B-spline knots}
#' \item{\code{len}}{curve length}
#' \item{\code{len_motifs}}{motif length}
#' \item{\code{N}}{number of curves}
#' \item{\code{nmotifs}}{number of motifs}
#' \item{\code{norder}}{order of B-splines}
#' \item{\code{sd_noise}}{standard deviation of Gaussian noise}
#' \item{\code{curves}}{list with the following elements:}
#' \item{\code{fd_curves}}{list of 18 functional data objects (class fd from
#' package fda) with the 18 curves}
#' \item{\code{motifs_in_curves}}{list of 18 objects, one for each curve,
#' reporting the ids of the motifs (id_motifs) present in the curve and their
#' starting position (pos_motifs) with respect to the B-spline basis (for
#' example a starting position equal to 3 means that the motif starts at the
#' third)}
#' \item{\code{coeff_motifs}}{list of 2 vectors with the B-spline
#' coefficients defining the two motifs}
#' \item{\code{freq_motifs}}{list of 18 vectors with the frequency (number of
#' occurrences) of the two motifs in each of the 18 curves}
#' }
#' @name foodinfl_loc
#' @usage data(foodinfl_loc)
NULL