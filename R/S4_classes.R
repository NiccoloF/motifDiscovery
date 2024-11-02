
setClassUnion("charOrNum", c("character", "numeric"))

#' @title motifSimulationS4Class
#' 
#' @description 
#' The `motifSimulation` class is an S4 class designed for simulating functional data with embedded motifs. This class is essential for modeling various aspects of the data, such as the number of curves, motif details, curve characteristics, and knot-based spline definitions. It allows users to generate realistic synthetic functional data to test and benchmark motif discovery algorithms like ProbKMA and funBIalign.
#' 
#' @slot N A numeric value representing the number of curves to be simulated. 
#' This parameter controls the number of functional curves generated in the simulation.
#' 
#' @slot mot_details A list containing details of the motifs to be embedded within the curves. 
#' Each motif detail specifies attributes such as motif shape, location, and frequency of occurrence across the simulated curves.
#' 
#' @slot motifs_in_curves A list specifying which curves contain motifs and the positions where they are embedded. 
#' This slot allows precise control over the placement of motifs in the generated curves, enabling flexible motif-to-curve assignments.
#' 
#' @slot distribution A character string or numeric value representing the distribution of the weights for the motifs. 
#' This slot defines the distribution used to generate the weight coefficients for the motifs, influencing their amplitude in the simulated data. 
#' Accepted values include any distribution supported in R or a custom vector provided by the user.
#' 
#' @slot dist_knots A numeric value indicating the distance between knots in the spline representation of the curves. 
#' This parameter is crucial for defining the smoothness of the generated curves and the precision of the spline-based motif embedding process.
#' 
#' @slot len A numeric value representing the length of the generated curves. 
#' This parameter determines the number of time points or the granularity of the data, allowing for high-resolution simulations.
#' 
#' @slot norder A numeric value specifying the order of the B-spline used in the functional data representation. 
#' Higher-order splines provide smoother representations of the curves, while lower-order splines offer more flexible curve shapes.
#' 
#' @slot coeff_min A numeric value specifying the minimum coefficient value for the spline-based curve generation. 
#' This controls the lower bound of the weight coefficients applied to the basis functions, influencing the range of curve amplitudes.
#' 
#' @slot coeff_max A numeric value specifying the maximum coefficient value for the spline-based curve generation. 
#' This controls the upper bound of the weight coefficients applied to the basis functions, setting the maximum amplitude for the curves.
#' 
#' @slot min_dist_motifs A numeric value indicating the minimum distance between motifs in the simulated curves. 
#' This ensures that motifs are not placed too closely, preserving their distinctiveness and reducing overlap during the simulation.
#' 
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