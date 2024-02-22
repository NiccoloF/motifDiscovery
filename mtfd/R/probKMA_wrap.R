#' @title probKMA_wrap
#'
#' @description R wrap for function probKMA used to parallelize find_candidate_motifs
#' @param Y0 list of N vectors, for univariate curves y_i(x), or list of N matrices with d columns,
#' for d-dimensional curves y_i(x), with the evaluation of curves (all curves should be evaluated on
#' a uniform grid). When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
#' @param Y1 list of N vectors, for univariate derivative curves y'_i(x), or list of N matrices with
#' d columns, for d-dimensional derivatibe curves y'_i(x), with the evaluation of the curves derivatives
#' (all curves should be evaluated on a uniform grid). When y'_j(x)=NA in the dimension j, then y'_j(x)=NA
#' in ALL dimensions. Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
#' @param standardize if TRUE, each dimension is standardized (Z score on all the regions together).
#' @param K number of motifs.
#' @param c minimum motif lengths. Can be an integer (or a vector of K integers).
#' @param c_max maximum motif lengths. Can be an integer (or a vector of K integers).
#' @param P0 initial membership matrix, with N row and K column (if NULL, a random P0 is choosen).
#' @param S0 initial shift warping matrix, with N row and K column (if NULL, a random S0 is choosen).
#' @param diss dissimilarity. Possible choices are 'd0_L2', 'd1_L2', 'd0_d1_L2'.
#' @param alpha when diss='d0_d1_L2', weight coefficient between d0_L2 and d1_L2 (alpha=0 means d0_L2,
#' alpha=1 means d1_L2).
#' @param w vector of weights for the dissimilarity index in the different dimensions (w>0).
#' @param m m>1 weighting exponent in least-squares functional.
#' @param iter_max the maximum number of iterations allowed.
#' @param stop_criterion criterion to stop iterate, based on the Bhattacharyya distance
#' between memberships in subsequent iterations. Possible choices are: 'max' for the
#' maximum of distances in the different motifs; 'mean' for the average of distances in
#' the different motifs; 'quantile' for the quantile of distances in the different motifs
#' (in this case, quantile must be provided).
#' @param quantile probability in (0,1) to be used if stop.criterion='quantile'.
#' @param tol method tolerance (method stops if the stop criterion <tol).
#' @param iter4elong motifs elongation is performed every iter4elong iterations (if iter4elong>iter.max,
#' no elongation is done).
#' @param tol4elong tolerance on the Bhattacharyya distance (with the choice made in stop.criterion) for
#' performing motifs elongation.
#' @param max_elong maximum elongation allowed in a single iteration, as percentage of the motif length.
#' @param trials_elong number of (equispaced) elongation trials at each side of the motif in a single
#' iteration.
#' @param deltaJk_elong maximum relative objective function increasing allowed in each motif elongation
#' (for gaps and each side).
#' @param max_gap maximum gap allowed in each alignment (percentage of motif length).
#' @param iter4clean motif cleaning is performed every iter4clean iterations (if iter4clean>iter_max,
#' no cleaning is done).
#' @param tol4clean tolerance on the Bhattacharyya distance (with the choice made in stop_criterion) for
#' performing motif cleaning.
#' @param quantile4clean dissimilarity quantile to be used in motif cleaning.
#' @param return_options if TRUE, the options K,c,diss,w,m are returned by the function.
#' @param return_init if TRUE, P0 and S0 are returned by the function.
#' @param worker_number number of CPU cores to be used for parallelization (default number of CPU cores -1).
#' If worker_number=1, the function is run sequentially.
#' @return A list containing some initialization and input options of the function: P0, S0 (if return_init = TRUE), Y0, Y1, standardize, K, c, c_max, diss, alpha, w, m, iter4elong, tol4elong, max_elong, trials_elong, deltaJk_elong, max_gap, iter4clean and tol4clean (if return_options = TRUE);
#' @return \item{P}{ membership probability matrix}
#' @return \item{P_clean}{ membership probability matrix dichotomized according to the quantile of order 1/K}
#' @return \item{S}{ shift warping matrix}
#' @return \item{S_clean}{ shift warping matrix after cleaning motifs}
#' @return \item{D}{ dissimilarity matrix}
#' @return \item{D_clean}{ dissimilarity matrix after cleaning motifs}
#' @return \item{iter}{ iterations number}
#' @return \item{j_iter}{ minimum objective function}
#' @return \item{BC_dist}{ Bhattacharyya distance}
#' @return \item{BC_dist_iter}{ Bhattacharyya distance from the last iteration}
#' @return \item{V0}{ list of the candidates motifs}
#' @return \item{V0_clean}{ list of the candidates motifs after cleaning}
#' @return \item{V1}{ list of derived from candidates motifs}
#' @return \item{V1_clean}{ list of derived from candidates motifs after cleaning}
#' @author Riccardo Lazzarini Niccolo' Feresini
#' @export
probKMA_wrap <- function(Y0 = NULL,Y1 = NULL,P0 = matrix(),S0 = matrix(),
                         standardize= FALSE,c_max = Inf,iter_max = 1000,
                         iter4elong = 10,trials_elong = 10,return_options = TRUE,
                         alpha = 0,max_gap = 0.2,quantile = 0.25, stopCriterion = 'max', 
                         tol = 1e-8, tol4elong = 1e-3, max_elong = 0.5, deltaJK_elong = 0.05, 
                         iter4clean = 50, tol4clean = 1e-4,m = 2,w = 1, seed = 1, 
                         K = 2, c = 40, quantile4clean = 1/K, exe_print = FALSE,
                         set_seed = FALSE,n_threads = 7,diss = 'd0_2'){
  
  params = list(standardize=standardize,c_max = c_max,iter_max = iter_max,
                iter4elong = iter4elong,trials_elong = trials_elong,
                return_options = return_options, alpha = alpha,
                max_gap = max_gap,quantile = quantile, 
                stopCriterion = stopCriterion, tol = tol, 
                tol4elong = tol4elong, max_elong = max_elong, 
                deltaJK_elong = deltaJK_elong, iter4clean = iter4clean, 
                tol4clean = tol4clean,quantile4clean = quantile4clean, 
                m = m, w = w, seed = seed, K = K, c = c, exe_print = exe_print,
                set_seed = set_seed,n_threads = n_threads) 
  
  checked_data <- initialChecks(Y0,Y1,P0,S0,params,diss,seed)
  
  params <- checked_data$Parameters
  
  data <- checked_data$FuncData
  
  if ( alpha == 0 ||  alpha == 1)
  {
    string  = "L2"
  } 
  else 
  {
    string = "H1"
  }
  
  prok = new(ProbKMA,data$Y,params,data$P0,data$S0,string)
  
  rm(params)
  
  probKMA_results_1 = list(Y0 = data$Y$Y0,Y1 = data$Y$Y1,
                           diss = diss, w = w, alpha = alpha)
  
  rm(data)
  
  probKMA_results_2 = prok$probKMA_run() 
  
  rm(prok)
  
  return(c(probKMA_results_1,probKMA_results_2))
}
