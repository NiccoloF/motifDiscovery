#' @title MyS4Class
#' @description A simple S4 class to demonstrate S4 objects in R packages.
#' @slot name A character string representing the name.
#' @slot value A numeric value.
#' @export
setClassUnion("charOrNum", c("character", "numeric"))

#' @export
setClass("motifSimulation",
         slots = list(
           N = "numeric",            # Number of curves
           mot_details = "list",     # Details of motifs 
           motifs_in_curves = "list",# Details of curves
           distribution = "charOrNum",  # distribution of the weights
           dist_knots = "numeric",   # Distance between knots 
           len = "numeric",          # Length of the curve or curves
           norder = "numeric",       # Order of spline 
           coeff_min = "numeric",    # Minimum coefficient value
           coeff_max = "numeric",    # Maximum coefficient value
           min_dist_motifs = "numeric" # Minimum distance between motifs
         ))