#' @docType class
#' @title ProbKMA Class
#'
#' @description
#' The `ProbKMA` class is an R wrapper for the C++ implementation of the Probabilistic K-means Algorithm (ProbKMA) with local alignment.
#' This class facilitates local clustering of functional data and functional motif discovery, as proposed in the paper
#' `Probabilistic K-means with local alignment for clustering and motif discovery in functional data`,
#' authored by Marzia A. Cremona and Francesca Chiaromonte.
#'
#' Create a `ProbKMA` object from the ProbKMA C++ class.
#'
#' @return
#' A `ProbKMA` object from the C++ ProbKMA class.
#'
#' @param Y A list containing functional data and possibly derivatives.
#' @param Parameters An instance of the Parameters class, containing algorithm settings.
#' @param P0 A matrix representing the initial membership probabilities.
#' @param S0 A matrix representing the initial shift warping parameters.
#' @param diss A character string specifying the dissimilarity measure. Possible choices are:
#'   - `'d0_L2'`
#'   - `'d1_L2'`
#'   - `'d0_d1_L2'`
#'
#' @examples
#' \dontrun{
#' # Seed for random initialization of P0 and S0
#' seed <- 1
#'
#' # Set type of distance
#' diss <- 'd0_d1_L2'  # options: 'd0_L2', 'd1_L2', 'd0_d1_L2'
#'
#' # Null matrix for random initialization
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
#' a <- ProbKMAcpp::initialChecks(simulated200$Y0, simulated200$Y1, P0, S0, params, diss, seed)
#'
#' # Get data and parameters
#' params <- a$Parameters
#' data <- a$FuncData
#'
#' # Create an object of the class ProbKMA
#' prok <- new(ProbKMAcpp::ProbKMA, data$Y, data$V, params, data$P0, data$S0, "H1")
#'
#' # Run ProbKMA algorithm
#' output <- prok$probKMA_run()
#' }
#'
#' ##################
#' ## Constructor
#' prok <- new(ProbKMA, data$Y, data$V, params, data$P0, data$S0, "H1")  # or 'L2'
#'
#' ##################
#' ## Getters
#' prok$get_parameters()  # returns a list of parameters
#' prok$get_motifs()      # returns a list containing the motifs found
#'
#' ##################
#' ## Setters
#' prok$set_P0(P)         # set membership matrix
#' prok$set_S0(S)         # set shift warping matrix
#' prok$set_parameters(param)  # set parameters field by passing a list of parameters
#'
#' ##################
#' ## Initialize Motifs
#' prok$reinit_motifs(c, d)  # reinitialize (empty) K motifs with dimension c_k x d
#'
#' ##################
#' ## Run ProbKMA algo
#' prok$probKMA_run()  # run the algorithm
#'
#' @author Niccolò Feresini and Riccardo Lazzarini
#' @name ProbKMA
#' @export ProbKMA
#'
#' # Load the Rcpp module exposed with RCPP_MODULE( ... ) macro.
loadModule(module = "ProbKMAModule", TRUE)