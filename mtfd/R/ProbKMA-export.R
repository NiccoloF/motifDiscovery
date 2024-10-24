#' @docType class
#' @title ProbKMA Class
#' @description 
#' The `ProbKMA` class is an R wrapper for the C++ implementation of the 
#' Probabilistic K-means Algorithm (ProbKMA) with local alignment. This class 
#' facilitates local clustering of functional data and functional motif discovery, 
#' as proposed in the paper 
#' `Probabilistic K-means with local alignment for clustering and motif discovery in functional data`, 
#' authored by Marzia A. Cremona and Francesca Chiaromonte.
#' 
#' @section Constructor:
#' Create a `ProbKMA` object using the following constructor:
#' 
#' \preformatted{
#' prok <- new(ProbKMA, data$Y, data$V, params, data$P0, data$S0, "H1")
#' }
#' 
#' @return 
#' A `ProbKMA` object from the C++ ProbKMA class.
#' 
#' @section Parameters:
#' \describe{
#'   \item{Y}{A list containing functional data and possibly derivatives.}
#'   \item{params}{An instance of the Parameters class, containing algorithm settings.}
#'   \item{P0}{A matrix representing the initial membership probabilities.}
#'   \item{S0}{A matrix representing the initial shift warping parameters.}
#'   \item{diss}{A character string specifying the dissimilarity measure. Possible choices are: 
#'     \itemize{
#'       \item `'d0_L2'`
#'       \item `'d1_L2'`
#'       \item `'d0_d1_L2'`
#'     }}
#' }
#' 
#' @section Usage:
#' You can access and modify the `ProbKMA` object with the following methods:
#' \describe{
#'   \item{Getters:}{
#'     \describe{
#'       \item{prok$get_parameters()}{Returns a list of parameters.}
#'       \item{prok$get_motifs()}{Returns a list containing the motifs found.}
#'     }
#'   }
#'   \item{Setters:}{
#'     \describe{
#'       \item{prok$set_P0(P)}{Sets the membership matrix.}
#'       \item{prok$set_S0(S)}{Sets the shift warping matrix.}
#'       \item{prok$set_parameters(param)}{Sets parameters field by passing a list of parameters.}
#'     }
#'   }
#'   \item{Initialize Motifs:}{
#'     \describe{
#'       \item{prok$reinit_motifs(c, d)}{Reinitializes (empty) K motifs with dimension c_k x d.}
#'     }
#'   }
#'   \item{Run ProbKMA algorithm:}{
#'     \describe{
#'       \item{prok$probKMA_run()}{Runs the algorithm.}
#'     }
#'   }
#' }
#' 
#' @examples
#' \dontrun{
#' # Example usage
#' # Seed for random initialization of P0 and S0
#' seed <- 1
#' 
#' # Set type of distance
#' diss <- 'd0_d1_L2'  # options: 'd0_L2', 'd1_L2', 'd0_d1_L2'
#' 
#' # Null matrices for random initialization
#' P0 <- matrix()
#' S0 <- matrix()
#' 
#' # Define parameters for the algorithm
#' params <- list(standardize = TRUE, K = 2, c = 61, c_max = 71,
#'                iter_max = 1000, quantile = 0.25,
#'                stopCriterion = 'max', tol = 1e-8,
#'                iter4elong = 1, tol4elong = 1e-3, max_elong = 0.5,
#'                trials_elong = 201, deltaJK_elong = 0.05, max_gap = 0,
#'                iter4clean = 50, tol4clean = 1e-4,
#'                quantile4clean = 1/2, return_options = TRUE,
#'                m = 2, w = 1, alpha = 0.5, seed = seed, exe_print = TRUE,
#'                set_seed = TRUE)
#' 
#' # Check input data
#' a <- initialChecks(simulated200$Y0, simulated200$Y1, P0, S0, params, diss, seed)
#' 
#' # Get data and parameters
#' params <- a$Parameters
#' data <- a$FuncData
#' 
#' # Create an object of the class ProbKMA
#' prok <- new(ProbKMA, data$Y, data$V, params, data$P0, data$S0, "H1")
#' 
#' # Run ProbKMA algorithm
#' output <- prok$probKMA_run()
#' }
#' 
#' @author Niccolò Feresini and Riccardo Lazzarini
#' @name ProbKMA
#' @export ProbKMA
loadModule(module = "ProbKMAModule", TRUE)
