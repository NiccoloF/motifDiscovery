#' @title Wrapper for the Probabilistic K-means Algorithm (ProbKMA)
#'
#' @description This function serves as a wrapper for the Probabilistic K-means Algorithm (ProbKMA) to cluster functional data. It handles preprocessing, parameter setup, and execution of the core algorithm, returning the results along with silhouette analysis to assess the clustering quality.
#'
#' @param Y0 A matrix of functional data for the first set of observations.
#' @param Y1 A matrix of functional data for the second set of observations.
#' @param P0 A matrix representing the initial membership probabilities.
#' @param S0 A matrix representing the initial shift parameters.
#' @param standardize A logical value indicating whether to standardize the data. Default is `FALSE`.
#' @param c_max Maximum number of motifs to extract. Default is `Inf`.
#' @param iter_max Maximum number of iterations for the algorithm. Default is 1000.
#' @param iter4elong Number of iterations for elongation. Default is 10.
#' @param trials_elong Number of trials for elongation. Default is 10.
#' @param return_options A logical value indicating whether to return additional options. Default is `TRUE`.
#' @param alpha A numeric value representing the weighting parameter. Default is 0.
#' @param max_gap Maximum allowable gap between motifs. Default is 0.2.
#' @param quantile Quantile to be used for cleaning. Default is 0.25.
#' @param stopCriterion Stopping criterion for the algorithm, can be 'max' or other specified values. Default is 'max'.
#' @param tol Tolerance for convergence. Default is 1e-8.
#' @param tol4elong Tolerance for elongation iterations. Default is 1e-3.
#' @param max_elong Maximum elongation allowed. Default is 0.5.
#' @param deltaJK_elong Increment for the elongation. Default is 0.05.
#' @param iter4clean Number of iterations for the cleaning process. Default is 50.
#' @param tol4clean Tolerance for the cleaning process. Default is 1e-4.
#' @param m Parameter controlling the clustering behavior. Default is 2.
#' @param w Weighting parameter for the dissimilarity measure. Default is 1.
#' @param seed Random seed for reproducibility. Default is 1.
#' @param K Number of motifs to extract. Default is 2.
#' @param c Minimum motif length. Default is 40.
#' @param quantile4clean Quantile used for the cleaning process. Default is 1/K.
#' @param exe_print A logical value indicating whether to print execution details. Default is `FALSE`.
#' @param set_seed A logical value indicating whether to set the random seed. Default is `FALSE`.
#' @param diss Dissimilarity measure to be used. Default is 'd0_2'.
#' @param transformed A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.
#' @param V_init Initial values for the motifs. Default is `NULL`.
#' @param align A logical value indicating whether to align the curves. Default is `TRUE`.
#' @param n_threads Number of threads to use for parallel processing. Default is 1.
#'
#' @return A list containing:
#' \item{probKMA_results}{A list of results from the ProbKMA algorithm, including processed functional data and model parameters.}
#' \item{silhouette_results}{Results from silhouette analysis, indicating the quality of the clustering.}
#'
#' @export
probKMA_wrap <- function(Y0 = NULL,Y1 = NULL,P0 = matrix(),S0 = matrix(),
                         standardize= FALSE,c_max = Inf,iter_max = 1000,
                         iter4elong = 10,trials_elong = 10,return_options = TRUE,
                         alpha = 0,max_gap = 0.2,quantile = 0.25, stopCriterion = 'max', 
                         tol = 1e-8, tol4elong = 1e-3, max_elong = 0.5, deltaJK_elong = 0.05, 
                         iter4clean = 50, tol4clean = 1e-4,m = 2,w = 1, seed = 1, 
                         K = 2, c = 40, quantile4clean = 1/K, exe_print = FALSE,
                         set_seed = FALSE,diss = 'd0_2', 
                         transformed = FALSE, V_init = NULL,align = TRUE,n_threads = 1){
  params = list(standardize=standardize,c_max = c_max,iter_max = iter_max,
                iter4elong = iter4elong,trials_elong = trials_elong,
                return_options = return_options, alpha = alpha,
                max_gap = max_gap,quantile = quantile, 
                stopCriterion = stopCriterion, tol = tol, 
                tol4elong = tol4elong, max_elong = max_elong, 
                deltaJK_elong = deltaJK_elong, iter4clean = iter4clean, 
                tol4clean = tol4clean,quantile4clean = quantile4clean, 
                m = m, w = w, seed = seed, K = K, c = c, exe_print = exe_print,
                set_seed = set_seed, 
                transformed = transformed, V_init = V_init,n_threads = n_threads) 
  checked_data <- initialChecks(Y0,Y1,P0,S0,params,diss,V_init)
  params <- checked_data$Parameters
  
  data <- checked_data$FuncData
  string_diss <- ifelse(alpha == 0 || alpha == 1,"L2","H1")
  prok <- NULL
  if (!is.null(data$V_init)) {
    prok <- new(ProbKMA, data$Y, params, data$P0, data$S0, string_diss, data$V_init)
  } else {
    prok <- new(ProbKMA, data$Y, params, data$P0, data$S0, string_diss)
  }
  probKMA_results_1 = list(Y0 = data$Y$Y0,
                           Y1 = data$Y$Y1,
                           diss = diss, w = w, 
                           alpha = alpha, 
                           V_init = data$V_init, 
                           transformed = transformed)
  
  probKMA_results_2 = prok$probKMA_run() 
  sil <- prok$compute_silhouette(align)
  
  rm(params)
  rm(data)
  rm(prok)
  return(list(probKMA_results = c(probKMA_results_1,probKMA_results_2), silhouette_results = sil))

}
