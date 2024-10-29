#' @title Functional Motif Discovery
#'
#' @description
#' The `discoverMotifs` function facilitates the discovery of recurring patterns, or motifs, within functional data by employing two sophisticated algorithms: 
#' \code{ProbKMA} (Probabilistic K-means with Local Alignment) and \code{funBIalign}. 
#' These algorithms are designed to identify and cluster functional motifs across multiple curves, leveraging advanced clustering and alignment techniques to handle complex data structures.
#'
#' \code{ProbKMA} integrates probabilistic clustering with local alignment strategies, enabling the detection of motifs that exhibit variability in both shape and position across different curves. 
#' This method is particularly adept at handling noisy data and motifs that may appear at varying scales or locations within the curves.
#'
#' On the other hand, \code{funBIalign} utilizes hierarchical clustering based on mean squared residue scores to uncover motifs. 
#' This approach effectively captures the additive nature of functional motifs, considering both portion-specific adjustments and time-varying components to accurately identify recurring patterns.
#'
#' By providing a flexible interface that accommodates different clustering paradigms, `discoverMotifs` empowers users to perform robust motif discovery tailored to their specific data characteristics and analytical requirements. 
#' Whether opting for the probabilistic and alignment-focused \code{ProbKMA} or the hierarchical and residue-based \code{funBIalign}, users can leverage these methods to extract meaningful and interpretable motifs from their functional datasets.
#'
#' @details
#' The `discoverMotifs` function dynamically switches between two advanced motif discovery algorithms based on the user's specification. 
#' Each algorithm employs distinct strategies to identify and cluster motifs within functional data, offering flexibility and adaptability to various analytical scenarios.
#'
#' @section Theoretical Background for ProbKMA:
#' 
#' \code{ProbKMA} is inspired by methodologies prevalent in bioinformatics, particularly those involving local alignment techniques extended from high-similarity seeds. 
#' This algorithm combines fuzzy clustering approaches with local alignment strategies to effectively minimize a generalized least squares functional. 
#' The minimization process can incorporate both the levels and derivatives of the curves through a Sobolev-based distance metric, enhancing the algorithm's sensitivity to both shape and rate changes in the data.
#' 
#' Throughout its iterative process, \code{ProbKMA} refines motif centers, membership probabilities, and alignment shifts, making it highly effective for capturing complex motif structures and motifs distributed across multiple curves. 
#' This ensures that the discovered motifs are both representative and robust against variations and noise within the functional data.
#'
#' @section Theoretical Background for funBIalign:
#' 
#' \code{funBIalign} models functional motifs as an additive combination of motif means, portion-specific adjustments, and time-varying components. 
#' The algorithm constructs a hierarchical dendrogram utilizing the generalized mean squared residue score (fMSR) to identify candidate motifs across curves.
#' 
#' A critical aspect of \code{funBIalign} is its post-processing step, which filters out redundant motifs and refines the final selection to ensure that only the most significant and representative motifs are retained. 
#' This hierarchical approach allows for a nuanced identification of motifs, capturing both broad and subtle patterns within the data.
#'
#' @section Common Parameters:
#' The following parameters are common to both \code{ProbKMA} and \code{funBIalign} algorithms:
#' \describe{
#'   \item{\code{Y0}}{A list containing \code{N} vectors (for univariate curves) or \code{N} matrices (for multivariate curves) representing the functional data. Each curve is evaluated on a uniform grid, ensuring consistency across the dataset.}
#'   \item{\code{method}}{A character string specifying the motif discovery algorithm to use. Acceptable values are \code{"ProbKMA"} for Probabilistic K-means with Local Alignment and \code{"funBIalign"} for Functional Bi-directional Alignment.}
#'   \item{\code{stopCriterion}}{A character string indicating the convergence criterion for the selected algorithm. For \code{ProbKMA}, options include \code{"max"}, \code{"mean"}, or \code{"quantile"} based on the Bhattacharyya distance between memberships in successive iterations. For \code{funBIalign}, options are \code{"fMRS"} (functional Mean Squared Residue) or \code{"Variance"} to guide the ranking of motifs.}
#'   \item{\code{name}}{A character string specifying the name of the output directory where results will be saved. This facilitates organized storage and easy retrieval of analysis results.}
#'   \item{\code{plot}}{A logical value indicating whether to generate and save plots of the discovered motifs and clustering results. When set to \code{TRUE}, visualizations are produced to aid in the qualitative assessment of the motif discovery process.}
#'   \item{\code{worker_number}}{An integer specifying the number of CPU cores to utilize for parallel computations. By default, the function uses the total number of available cores minus one, optimizing computational efficiency without overloading the system.}
#' }
#'
#' @section ProbKMA Options:
#' The following parameters are specific to the \code{ProbKMA} algorithm:
#' \describe{
#'   \item{\code{K}}{An integer or vector specifying the number of motifs to be discovered. It can be a single integer for uniform motif discovery or a vector for specifying different numbers of motifs.}
#'   \item{\code{c}}{An integer or vector indicating the minimum motif lengths. This ensures that each discovered motif meets a specified minimum length requirement, maintaining the integrity of motif structures.}
#'   \item{\code{c_max}}{An integer or vector specifying the maximum motif lengths, allowing control over the upper bounds of motif sizes to prevent excessively long motifs.}
#'   \item{\code{diss}}{A character string defining the dissimilarity measure to use. Possible values include \code{"d0_L2"}, \code{"d1_L2"}, and \code{"d0_d1_L2"}, which determine how the algorithm quantifies differences between motifs based on level and derivative information.}
#'   \item{\code{alpha}}{A numeric value between 0 and 1 that serves as a weight parameter between \code{d0_L2} and \code{d1_L2} when using \code{d0_d1_L2}. An \code{alpha} of 0 emphasizes \code{d0_L2}, while an \code{alpha} of 1 emphasizes \code{d1_L2}, allowing for balanced consideration of both metrics.}
#'   \item{\code{w}}{A numeric vector specifying the weight for the dissimilarity index across different dimensions. All values must be positive, enabling the algorithm to prioritize certain dimensions over others based on their relative importance.}
#'   \item{\code{m}}{A numeric value greater than 1 that acts as the weighting exponent in the least-squares functional method. This parameter influences the sensitivity of the algorithm to differences in motif alignment and membership probabilities.}
#'   \item{\code{iter_max}}{An integer specifying the maximum number of iterations allowed for the algorithm to converge. This prevents excessive computation time by limiting the number of optimization steps.}
#'   \item{\code{quantile}}{A numeric value representing the quantile probability used when \code{stopCriterion} is set to \code{"quantile"}. This determines the threshold for convergence based on the distribution of Bhattacharyya distances.}
#'   \item{\code{tol}}{A numeric value specifying the tolerance level for convergence. The algorithm stops iterating if the change in the stop criterion falls below this threshold, ensuring precise and stable convergence.}
#'   \item{\code{iter4elong}}{An integer indicating the number of iterations after which motif elongation is performed. If set to a value greater than \code{iter_max}, no elongation is performed. Motif elongation allows the algorithm to extend motifs to better fit the data.}
#'   \item{\code{tol4elong}}{A numeric value defining the tolerance on the Bhattacharyya distance for motif elongation. This parameter controls how much the objective function can increase during elongation, ensuring that motif extensions do not degrade the overall fit.}
#'   \item{\code{max_elong}}{A numeric value representing the maximum elongation allowed in a single iteration, expressed as a percentage of the motif length. This prevents excessive extension of motifs in any single step.}
#'   \item{\code{trials_elong}}{An integer specifying the number of elongation trials (equispaced) on each side of the motif in a single iteration. Multiple trials enhance the robustness of motif elongation by exploring various extension possibilities.}
#'   \item{\code{deltaJK_elong}}{A numeric value indicating the maximum relative increase in the objective function permitted during motif elongation. This ensures that elongation steps contribute positively to the motif fitting process.}
#'   \item{\code{max_gap}}{A numeric value defining the maximum gap allowed in each alignment as a percentage of the motif length. This parameter controls the allowable discontinuity between aligned motifs, maintaining coherence in motif placement.}
#'   \item{\code{iter4clean}}{An integer specifying the number of iterations after which motif cleaning is performed. If set to a value greater than \code{iter_max}, no cleaning is performed. Motif cleaning removes redundant or poorly fitting motifs to refine the final motif set.}
#'   \item{\code{tol4clean}}{A numeric value representing the tolerance on the Bhattacharyya distance for motif cleaning. This parameter determines the threshold for identifying and removing redundant motifs during the cleaning process.}
#'   \item{\code{quantile4clean}}{A numeric value specifying the dissimilarity quantile used for motif cleaning. This quantile determines which motifs are considered sufficiently dissimilar to be retained in the final set.}
#'   \item{\code{return_options}}{A logical value indicating whether to return the options passed to the \code{ProbKMA} method. When set to \code{TRUE}, users receive detailed information about the algorithm's configuration, facilitating transparency and reproducibility.}
#'   \item{\code{Y1}}{A list of derivative curves used if the dissimilarity measure \code{"d0_d1_L2"} is selected. These derivatives enhance the algorithm's ability to capture both shape and rate changes in the functional data.}
#'   \item{\code{P0}}{An initial membership matrix (N x K), where N is the number of curves and K is the number of clusters. If set to \code{NULL}, a random matrix is generated, initiating the probabilistic clustering process.}
#'   \item{\code{S0}}{An initial shift warping matrix (N x K). If set to \code{NULL}, a random matrix is generated to initialize the alignment process, allowing motifs to adapt to variations in the data.}
#'   \item{\code{n_subcurves}}{An integer specifying the number of splitting subcurves used when the number of curves is equal to one. This parameter allows the algorithm to handle single-curve datasets by dividing them into manageable segments for motif discovery.}
#'   \item{\code{sil_threshold}}{A numeric value representing the threshold to filter candidate motifs based on their silhouette scores. This ensures that only motifs with sufficient clustering quality are retained in the final results.}
#'   \item{\code{set_seed}}{A logical value indicating whether to set a random seed for reproducibility. When set to \code{TRUE}, the function initializes the random number generator to ensure consistent results across multiple runs.}
#'   \item{\code{seed}}{An integer specifying the random seed used for initialization when \code{set_seed} is \code{TRUE}. This parameter guarantees reproducibility of the clustering and alignment processes.}
#'   \item{\code{exe_print}}{A logical value determining whether to print execution details for each iteration. When set to \code{TRUE}, users receive real-time feedback on the algorithm's progress, aiding in monitoring and debugging.}
#'   \item{\code{V_init}}{A list of motif sets provided as specific initializations for clustering rather than using random initializations. The `V_init` parameter allows users to provide a set of motifs as starting points for the algorithm, instead of relying on random initialization. If `n_init` is specified as greater than the number of motifs given in `V_init`, the remaining initializations will be randomly generated. For example, if `n_init = 10` but only 5 motif sets are given in `V_init`, the algorithm will use these 5 initializations and generate an additional 5 randomly.}
#'   \item{\code{transformed}}{A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.}
#'   \item{\code{n_init_motif}}{The number of initial motif sets from `V_init` to be used directly as starting points in clustering. If `n_init_motif` is set to a value larger than the number of motifs provided in `V_init`, additional initializations will be generated randomly to meet the specified number. For example, if `n_init = 10` and `n_init_motif = 5` with only 3 motif sets in `V_init`, the algorithm will use these 3 sets and generate 7 additional random initializations.}
#' }
#'
#' @section funBIalign Options:
#' The following parameters are specific to the \code{funBIalign} algorithm:
#' \describe{
#'   \item{\code{portion_len}}{An integer specifying the length of curve portions to align. This parameter controls the granularity of alignment, allowing the algorithm to focus on specific segments of the curves for motif discovery.}
#'   \item{\code{min_card}}{An integer representing the minimum cardinality of motifs, i.e., the minimum number of motif occurrences required for a motif to be considered valid. This ensures that only motifs with sufficient representation across the dataset are retained.}
#'   \item{\code{cut_off}}{A double that specifies the number of top-ranked motifs to keep based on the ranking criteria, facilitating focused visualization of the most significant motifs.
#'                         In particular, all motifs that rank below the cut_off are retained.}
#' }
#'
#' @param Y0 A list containing N vectors (for univariate curves) or N matrices (for multivariate curves) representing the functional data.
#' @param method A character string specifying the motif discovery algorithm to use. Acceptable values are "ProbKMA" for Probabilistic K-means with Local Alignment and "funBIalign" for Functional Bi-directional Alignment.
#' @param stopCriterion A character string indicating the convergence criterion for the selected algorithm.
#' @param name A character string specifying the name of the output directory where results will be saved.
#' @param plot A logical value indicating whether to generate and save plots of the discovered motifs and clustering results.
#' @param probKMA_options A list of options specific to the ProbKMA algorithm.
#' @param funBIalign_options A list of options specific to the funBIalign algorithm.
#' @param worker_number An integer specifying the number of CPU cores to utilize for parallel computations.
#'
#' @return 
#' A list containing the discovered motifs and their corresponding statistics, tailored to the selected method:
#' \describe{
#'   \item{\code{motifs}}{A list of identified motifs, each containing the motif's representative curve, membership probabilities, and alignment information.}
#'   \item{\code{statistics}}{Detailed statistics for each motif, including measures such as silhouette scores, variance explained, and other relevant metrics that quantify the quality and significance of the discovered motifs.}
#'   \item{\code{parameters}}{The final parameters and configurations used during the motif discovery process, providing transparency and facilitating reproducibility of the results.}
#'   \item{\code{plots}}{If \code{plot = TRUE}, this component contains the generated plots visualizing the motifs and their distribution across the functional data.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Discover motifs using ProbKMA
#' 
#' # Define dissimilarity measure and weight parameter
#' diss <- 'd0_d1_L2'
#' alpha <- 0.5
#' 
#' # Define number of motifs and their minimum lengths
#' K <- c(2, 3)
#' c <- c(61, 51)
#' n_init <- 10
#' 
#' # Load simulated data
#' data("simulated200")
#' 
#' # Perform motif discovery using ProbKMA
#' results <- mtfd::discoverMotifs(
#'   Y0 = simulated200$Y0,
#'   method = "ProbKMA",
#'   stopCriterion = "max",
#'   name = './results_ProbKMA_VectorData/',
#'   plot = TRUE,
#'   probKMA_options = list(
#'     Y1 = simulated200$Y1,
#'     K = K,
#'     c = c,
#'     n_init = n_init,
#'     diss = diss,
#'     alpha = alpha
#'   ),
#'   worker_number = NULL
#' )
#' 
#' # Modify silhouette threshold and re-run post-processing
#' results <- mtfd::discoverMotifs(
#'   Y0 = simulated200$Y0,
#'   method = "ProbKMA",
#'   stopCriterion = "max",
#'   name = './results_ProbKMA_VectorData/',
#'   plot = TRUE,
#'   probKMA_options = list(
#'     Y1 = simulated200$Y1,
#'     K = K,
#'     c = c,
#'     n_init = n_init,
#'     diss = diss,
#'     alpha = alpha,
#'     sil_threshold = 0.5
#'   ),
#'   worker_number = NULL
#' )
#' 
#' # Example 2: Discover motifs using funBIalign
#' results_funbialign <- mtfd::discoverMotifs(
#'   Y0 = simulated200$Y0,
#'   method = "funBIalign",
#'   stopCriterion = 'Variance',
#'   name = './results_FunBialign/',
#'   plot = TRUE,
#'   funBIalign_options = list(
#'     portion_len = 60,
#'     min_card = 3,
#'     cut_off = 1.0
#'   )
#' )
#' }
#'
#' @seealso 
#' \strong{ProbKMA}: 
#' \href{https://arxiv.org/pdf/1808.04773}{Probabilistic K-means with Local Alignment} \cr
#' \strong{funBIalign}: 
#' \href{https://arxiv.org/pdf/2306.04254}{Hierarchical Clustering with Mean Squared Residue Scores}.
#'
#' @export
discoverMotifs <- function(Y0,method,stopCriterion,name,plot,
                         probKMA_options = list(),
                         funBIalign_options = list(portion_len = NULL,min_card = NULL,cut_off=NULL),
                         worker_number = NULL){
  
  
  ### set parallel jobs #############################################################################
  core_number <- parallel::detectCores()
  # check worker number
  if(!is.null(worker_number)){
    if(!is.numeric(worker_number)){
      warning('Worker number not valid. Selecting default number.')
      worker_number=NULL
    }else{
      if((worker_number%%1!=0)|(worker_number<1)|(worker_number>core_number)){
        warning('Worker number not valid. Selecting default number.')
        worker_number=NULL
      }
    }
  } else {
    worker_number <- core_number-1 
  }
  rm(core_number)
  
  if(method == 'ProbKMA')
  {
    
    pb <- progress_bar$new(
      format = "Progress [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = 5,  
      width = 50
    )
    
    default_probKMA_options <- list(K=NULL,c=NULL,diss="d0_L2",alpha=0,
                                    Y1=NULL,P0=matrix(),S0=matrix(),
                                    n_init=10,return_options=TRUE,names_var='x(t)',
                                    standardize=FALSE,c_max=Inf,
                                    iter_max=1e3,iter4elong=100,
                                    trials_elong=201,max_gap=0.2,
                                    quantile=0.25,tol=1e-8,
                                    tol4elong=1e-3,max_elong=0.5,
                                    deltaJK_elong=0.05,iter4clean=50,
                                    tol4clean=1e-4,quantile4clean=0.5,
                                    m=2,w=1,transformed = FALSE,silhouette_align=FALSE,
                                    n_subcurves = 10,sil_threshold=0.9,
                                    seed = 1,exe_print = FALSE,set_seed = FALSE,
                                    V_init=NULL, n_init_motif = 0)
    probKMA_options <- modifyList(default_probKMA_options, probKMA_options)
    Y1 <- probKMA_options$Y1
    P0 <- probKMA_options$P0
    S0 <- probKMA_options$S0
    K <- probKMA_options$K
    c <- probKMA_options$c
    n_init <- probKMA_options$n_init
    diss <- probKMA_options$diss
    alpha <- probKMA_options$alpha
    return_options <- probKMA_options$return_options
    names_var <- probKMA_options$names_var
    seed <- probKMA_options$seed
    exe_print <- probKMA_options$exe_print
    set_seed <- probKMA_options$set_seed
    V_init <- probKMA_options$V_init
    n_init_motif <- probKMA_options$n_init_motif
    
    ### check input #############################################################################################
    # check name
    if(!is.character(name))
      stop('name should be a string.')
    if(grepl(' ',name)){
      warning('name contains spaces. Changing spacing in underscores.')
      name=gsub(" ","_",name,fixed=TRUE)
    }
    if (substr(name, nchar(name), nchar(name)) != "/") {
      name <- paste0(name, "/")
    }
    files = list.files(name)
    if('find_candidate_motifs_results.RData' %in% files) {
      # candidate motifs already present, load them
      load(paste0(name,'find_candidate_motifs_results.RData'))
    } else {
      
    ### check on n_init_motif (to be checked) ############
    if(n_init_motif > n_init)
    {
      warning('n_init_motif greater than n_init: setting n_init_motif = n_init')
      n_init_motif = n_init
    }
    if(n_init_motif%%1!=0)
      stop('Number of initializations n_init_motif should be an integer.')
    if(n_init_motif<0)
      stop('Number of initializations n_init_motif should be positive.')
    ######################################################
    
    # check required input
    if(missing(K))
      stop('K must be specified')
    if(missing(c))
      stop('c must be specified')
    # check K
    if(!is.vector(K))
      stop('K should be a vector.')
    # check c
    if(!is.vector(c))
      stop('c should be a vector.')
    # check n_init
    if(n_init%%1!=0)
      stop('Number of initializations n_init should be an integer.')
    if(n_init<1)
      stop('Number of initializations n_init should be at least 1.')
    
    
    #### new check V_init #########################################
    d = ncol(Y0[[1]]) 
    
    if(n_init_motif > 0)
    {
      for(i in 1:length(K))
      {
        for(j in 1:length(c))
        {
          if(length(V_init[[i]][[j]]) < n_init_motif) 
          {
            stop('More initial motifs needed')
          }
          else
          {
            if(length(V_init[[i]][[j]]) > n_init_motif)
            {
              warning('Initial motifs passed more than necessary. Resize V_init to the correct size!')
              V_init[[i]][[j]] <- V_init[[i]][[j]][1:n_init_motif]
            }
            for(e in 1:(n_init-n_init_motif))
            {
              V_init[[i]][[j]] <- append(V_init[[i]][[j]], list(NULL))
            }
            for(k in 1:n_init_motif)
            {
              if(length(V_init[[i]][[j]][[k]])!= K[i])
              {
                stop('Number of initial motifs differs from K')
              }
              else
              {
                for(h in 1:K[i])
                {
                  if(nrow(V_init[[i]][[j]][[k]][[h]]$v0)!= c[j] || ncol(V_init[[i]][[j]][[k]][[h]]$v0)!= d)
                  {
                    stop('Uncorrect dimensions of initial motifs provided')
                  }
                  else
                  {
                    if(probKMA_options$diss == "d0_d1_L2" || probKMA_options$diss == "d1_L2")
                    {
                      if(nrow(V_init[[i]][[j]][[k]][[h]]$v1)!= c[j] || ncol(V_init[[i]][[j]][[k]][[h]]$v0)!= d)
                      {
                        stop('Uncorrect dimensions of initial motifs provided')
                      }
                    }
                    V_init_bool = TRUE
                  }
                }
              }
            }
          }
        }
      }
    }
    
    ##################################
    # check Y0
    if(missing(Y0))
      stop('Y0 must be specified.')
    if(!is.null(Y0))
      if(!is.list(Y0))
        stop('Y0 should be a list of vectors or matrices.')
    if((FALSE %in% lapply(Y0,is.matrix))&&(FALSE %in% lapply(Y0,is.vector)))
      stop('Y0 should be a list of vectors or matrices.')
    N=length(Y0) # number of curves
    if(N!=1 && N<5)
      stop('More curves y_i(x) needed.')
    Y0=lapply(Y0,as.matrix)
    # set return_options=TRUE
    if(!is.null(probKMA_options$return_options)){
      if(probKMA_options$return_options==FALSE){
        warning('Setting return_option=TRUE')
        probKMA_options$return_options=TRUE
      }
    }
   
    arguments = list(Y0 = Y0, Y1 = Y1, P0 = P0, S0 = S0, 
                     standardize = probKMA_options$standardize, c_max = probKMA_options$c_max,
                     iter_max = probKMA_options$iter_max, iter4elong = probKMA_options$iter4elong, 
                     trials_elong = probKMA_options$trials_elong, return_options = return_options,
                     alpha = alpha, max_gap = probKMA_options$max_gap, 
                     quantile = probKMA_options$quantile, stopCriterion = stopCriterion, tol = probKMA_options$tol,
                     tol4elong =probKMA_options$tol4elong, 
                     max_elong = probKMA_options$max_elong, deltaJK_elong = probKMA_options$deltaJK_elong,
                     iter4clean = probKMA_options$iter4clean, 
                     tol4clean = probKMA_options$tol4clean, m = probKMA_options$m, w = probKMA_options$w,
                     seed = seed, K = NULL, c = NULL, 
                     quantile4clean = probKMA_options$quantile4clean,
                     exe_print = exe_print, set_seed = set_seed, 
                     diss = diss, transformed = probKMA_options$transformed, v_init = NULL,
                     align = probKMA_options$silhouette_align,
                     n_threads = worker_number)
      
    if(worker_number>1){
        cl_find=parallel::makeCluster(worker_number,timeout=60*60*24*30)
        parallel::clusterExport(cl_find,c('name','names_var',
                                          'probKMA_plot','probKMA_silhouette_plot',
                                          '.mapply_custom','.diss_d0_d1_L2','.domain',
                                          '.select_domain','.find_min_diss',
                                          'probKMA_wrap','arguments'),envir=environment()) 
        parallel::clusterCall(cl_find, function() {
          combinat::combn
        })
        on.exit(parallel::stopCluster(cl_find))
      }else{
        cl_find=NULL
      }
    pb$update(0.1)
    ### run probKMA ##########################################################################################
    i_c_K = expand.grid(seq_len(n_init),c,K)
    vector_seed = seq(1,length(i_c_K$Var1))
    if(n_init_motif > 0 && V_init_bool)
    {
      V_init_unlist=unlist(unlist(V_init,recursive=FALSE),recursive=FALSE)
      results=tryCatch({.mapply_custom(cl_find,function(K,c,i,v_init,small_seed){
        dir.create(paste0(name,"K",K,"_c",c),showWarnings=FALSE,recursive = TRUE)
        files=list.files(paste0(name,"K",K,"_c",c))
        if(paste0('random',i,'.RData') %in% files){
          load(paste0(name,"K",K,"_c",c,'/random',i,'.RData')) 
          message("K",K,"_c",c,'_random',i,' loaded')
          return(list(probKMA_results=probKMA_results,
                      time=time,silhouette=silhouette))
        }else{
          iter=iter_max=1
          message("K",K,"_c",c,'_random',i,' running')
          while(iter==iter_max){
            start=proc.time()
            if(arguments$set_seed)
            {
              arguments$seed = small_seed
              set.seed(small_seed)
            }
            small_seed = small_seed + 1
            arguments$K = K
            arguments$c = c
            arguments$quantile4clean = 1/K
            arguments$v_init = v_init
            results = do.call(probKMA_wrap,arguments)
            probKMA_results = results[[1]]
            silhouette_results = results[[2]]
            end=proc.time()
            time=end-start
            iter=probKMA_results$iter
            iter_max=arguments$iter_max
            if(iter==iter_max)
              warning('Maximum number of iteration reached. Re-starting.')
          }
          if(plot) {
          transformed=probKMA_results$transformed
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
          probKMA_plot(probKMA_results,ylab=names_var,plot = plot,cleaned=FALSE, transformed = arguments$transformed)
          dev.off()
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
          probKMA_plot(probKMA_results,ylab=names_var,sil_avg = silhouette_results[[4]],plot = plot,cleaned=TRUE, transformed = arguments$transformed) 
          dev.off()
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
          silhouette = probKMA_silhouette_plot(silhouette_results,K,plot = plot)
          dev.off()
          } else {
            silhouette = probKMA_silhouette_plot(silhouette_results,K,plot = FALSE)
          }
          save(probKMA_results,time,silhouette,
               file=paste0(name,"K",K,"_c",c,'/random',i,'.RData'))
          return(list(probKMA_results=probKMA_results,
                      time=time,silhouette=silhouette))
        }
      },i_c_K[,3],i_c_K[,2],i_c_K[,1],V_init_unlist,vector_seed,SIMPLIFY=FALSE)},error = function(e){
        stop(e$message)
      })
    }
    if(n_init_motif == 0 || !V_init_bool)
    {
      if(N == 1) {
        
        .split_curve_randomly <- function(curve, expected_motif_length, n_subcurves, randomness_factor = 0.3) {
          length_curve <- nrow(curve)  # Number of observations in the curve
          
          # Step 1: Randomly generate candidate split points
          candidate_splits <- sort(sample(seq(expected_motif_length, length_curve - expected_motif_length, 1), n_subcurves * 2))
          
          # Step 2: Adjust split points with randomness while respecting the minimum length requirement
          split_points <- c(1)  # Start point
          last_split <- 1
          
          for (i in candidate_splits) {
            # Check if we can add a split while respecting minimum length
            if (i - last_split >= expected_motif_length) {
              segment <- curve[last_split:i, ]
              valid_indices <- which(!is.na(as.matrix(segment)[, 1]))  # Indices of valid values
              
              if (length(valid_indices) > 0) {
                valid_segments <- split(valid_indices, cumsum(c(1, diff(valid_indices) != 1)))  # Split into contiguous valid segments
                valid_segment <- FALSE
                for (seg in valid_segments) {
                  if (length(seg) >= expected_motif_length) {
                    valid_segment <- TRUE
                    break
                  }
                }
                # Introduce randomness in adding split points
                if (valid_segment && runif(1) < randomness_factor) {
                  split_points <- c(split_points, i)
                  last_split <- i
                }
              }
            }
          }
          
          # Add the last segment if it meets the minimum length requirement
          if (length_curve - last_split >= expected_motif_length) {
            split_points <- c(split_points, length_curve)
          } else {
            split_points[length(split_points)] <- length_curve
          }
          
          split_points <- unique(split_points)
          
          # Step 3: Adjust the number of split points to match n_subcurves
          if (length(split_points) - 1 > n_subcurves) {
            # Too many splits, reduce by merging closest split points
            while (length(split_points) - 1 > n_subcurves) {
              gaps <- diff(split_points)
              min_gap_index <- which.min(gaps)
              split_points <- split_points[-(min_gap_index + 1)]  # Remove the smaller split
            }
          } else if (length(split_points) - 1 < n_subcurves) {
            # Too few splits, add new split points in the largest gaps between existing points
            additional_splits_needed <- n_subcurves - (length(split_points) - 1)
            while (additional_splits_needed > 0) {
              gaps <- diff(split_points)
              max_gap_index <- which.max(gaps)
              new_split <- round(mean(split_points[max_gap_index:(max_gap_index + 1)]))
              
              if (new_split - split_points[max_gap_index] >= expected_motif_length &&
                  split_points[max_gap_index + 1] - new_split >= expected_motif_length) {
                split_points <- sort(c(split_points, new_split))
                additional_splits_needed <- additional_splits_needed - 1
              } else {
                stop("Failed to find valid split points. Ensure each subcurve meets the minimum length requirement.")
              }
            }
          }
          
          # Step 4: Create subcurves based on valid split points
          subcurves <- list()
          for (i in 1:(length(split_points) - 1)) {
            subcurve <- curve[split_points[i]:split_points[i + 1], ]
            if (nrow(as.matrix(subcurve)) > 0) {
              subcurves[[paste0("c", i)]] <- subcurve
            }
          }
          
          return(list(subcurves = subcurves, split_points = split_points))
        }
        
        n_subcurves <- probKMA_options$n_subcurves 
        if(arguments$diss == 'd0_d1_L2') {
          split_results <- lapply(seq_len(n_init),function(null){.split_curve_randomly(arguments$Y0[[1]], max(c),n_subcurves, randomness_factor = 0.3)})
          Y0_subcurves <- lapply(split_results,function(splitted_curve){splitted_curve[[1]]})
          Y1_subcurves <- lapply(seq_len(n_init),function(j){
                            temp <- list()
                            for (i in seq_along(split_results$split_points)[-length(split_results$split_points)])
                              temp[[i]] <- split_results[[j]][[1]][pos[i]:pos[i + 1]]
                            return(temp)
                             })
        }
        else if(arguments$diss == 'd0_L2') {
          split_results <- lapply(seq_len(n_init),function(null){.split_curve_randomly(arguments$Y0[[1]], max(c),n_subcurves, randomness_factor = 0.3)})
          Y0_subcurves <- lapply(split_results,function(splitted_curve){splitted_curve[[1]]})
          Y1_subcurves <- vector("list",n_init)
        }
        else if(arguments$diss == 'd1_L2') {
          Y0_subcurves <- vector("list",n_init)
          split_results <- lapply(seq_len(n_init),function(null){.split_curve_randomly(arguments$Y1[[1]], max(c),n_subcurves, randomness_factor = 0.3)})
          Y1_subcurves <- lapply(split_results,function(splitted_curve){splitted_curve[[1]]})
        }
      results=tryCatch({.mapply_custom(cl_find, function(K,c,i,Y0_curve,Y1_curve,small_seed){ 
        dir.create(paste0(name,"K",K,"_c",c),showWarnings=TRUE,recursive = TRUE)
        files=list.files(paste0(name,"K",K,"_c",c))
        if(paste0('random',i,'.RData') %in% files){
          load(paste0(name,"K",K,"_c",c,'/random',i,'.RData'))
          message("K",K,"_c",c,'_random',i,' loaded')
          return(list(probKMA_results=probKMA_results,
                      time=time,silhouette=silhouette))
        }else{
          iter=iter_max=1
          message("K",K,"_c",c,'_random',i,' running')
          if(arguments$set_seed)
          {
            arguments$seed = small_seed
            set.seed(small_seed)
          }
          while(iter==iter_max){
            start=proc.time()
            small_seed = small_seed + 1
            arguments$Y0 <- Y0_curve
            arguments$Y1 <- Y1_curve
            arguments$K = K
            arguments$c = c
            arguments$quantile4clean = 1/K
            results = do.call(probKMA_wrap,arguments)
            probKMA_results = results[[1]]
            silhouette_results = results[[2]]
            end=proc.time()
            time=end-start
            iter=probKMA_results$iter
            iter_max=arguments$iter_max
            
            if(iter==iter_max)
              warning('Maximum number of iteration reached. Re-starting.')
          }
          if(plot) {
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
          probKMA_plot(probKMA_results,
                              plot = TRUE,
                              ylab=names_var,
                              sil_avg = silhouette_results[[4]],
                              cleaned=FALSE,
                              transformed = arguments$transformed)
          dev.off()
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
          probKMA_plot(probKMA_results,
                              plot = TRUE,
                              ylab=names_var,
                              sil_avg = silhouette_results[[4]],
                              cleaned=TRUE,
                              transformed = arguments$transformed)
          dev.off()

          pdf(paste0(name,"K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
          silhouette = probKMA_silhouette_plot(silhouette_results,K,plot = plot)
          dev.off()
          } else {
            silhouette = probKMA_silhouette_plot(silhouette_results,K,plot = FALSE)
          }

          probKMA_results["v_init"] <- NULL # delete v_init components
          save(probKMA_results,time,silhouette,
               file=paste0(name,"K",K,"_c",c,'/random',i,'.RData'))
           return(list(probKMA_results=probKMA_results,
                      time=time,silhouette=silhouette))
         }
      },i_c_K[,3],i_c_K[,2],i_c_K[,1],Y0_subcurves,Y1_subcurves,vector_seed,SIMPLIFY=FALSE)},error = function(e){
        stop(e$message)
      })
      } else {
      results=tryCatch({.mapply_custom(cl_find, function(K,c,i,small_seed){ 
        dir.create(paste0(name,"K",K,"_c",c),showWarnings=TRUE,recursive = TRUE)
        files=list.files(paste0(name,"K",K,"_c",c))
        if(paste0('random',i,'.RData') %in% files){
          load(paste0(name,"K",K,"_c",c,'/random',i,'.RData'))
          message("K",K,"_c",c,'_random',i,' loaded')
          return(list(probKMA_results=probKMA_results,
                      time=time,silhouette=silhouette))
        }else{
          iter=iter_max=1
          message("K",K,"_c",c,'_random',i,' running')
          if(arguments$set_seed)
          {
            arguments$seed = small_seed
            set.seed(small_seed)
          }
          while(iter==iter_max){
            start=proc.time()
            small_seed = small_seed + 1
            arguments$K = K
            arguments$c = c
            arguments$quantile4clean = 1/K
            results = do.call(probKMA_wrap,arguments)
            probKMA_results = results[[1]]
            silhouette_results = results[[2]]
            end=proc.time()
            time=end-start
            iter=probKMA_results$iter
            iter_max=arguments$iter_max
            
            if(iter==iter_max)
              warning('Maximum number of iteration reached. Re-starting.')
          }
          if(plot) {
            pdf(paste0(name,"K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
            probKMA_plot(probKMA_results,
                         plot = TRUE,
                         ylab=names_var,
                         sil_avg = silhouette_results[[4]],
                         cleaned=FALSE,
                         transformed = arguments$transformed)
            dev.off()
            pdf(paste0(name,"K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
            probKMA_plot(probKMA_results,
                         plot = TRUE,
                         ylab=names_var,
                         sil_avg = silhouette_results[[4]],
                         cleaned=TRUE,
                         transformed = arguments$transformed)
            dev.off()
            
            pdf(paste0(name,"K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
            silhouette = probKMA_silhouette_plot(silhouette_results,K,plot = plot)
            dev.off()
          } else {
            silhouette = probKMA_silhouette_plot(silhouette_results,K,plot = FALSE)
          }
          
          probKMA_results["v_init"] <- NULL # delete v_init components
          save(probKMA_results,time,silhouette,
               file=paste0(name,"K",K,"_c",c,'/random',i,'.RData'))
          return(list(probKMA_results=probKMA_results,
                      time=time,silhouette=silhouette))
        }
      },i_c_K[,3],i_c_K[,2],i_c_K[,1],vector_seed,SIMPLIFY=FALSE)},error = function(e){
        stop(e$message)
      })
    }
    }

    results=split(results,list(factor(i_c_K[,2],c),factor(i_c_K[,3],K)))
    results=split(results,rep(K,each=length(c)))
    
    ### plot silhouette average #################################################################################
    silhouette_average_sd=lapply(results,
                                 function(results){
                                   silhouette_average=lapply(results,
                                                             function(results){
                                                               silhouette_average=numeric(n_init)
                                                               silhouette_sd=numeric(n_init)
                                                               for(i in seq_len(n_init)){
                                                                 silhouette_average[i]=mean(results[[i]]$silhouette$silhouette_average)
                                                                 silhouette_sd[i]=sd(results[[i]]$silhouette$silhouette_average)
                                                               }
                                                               return(cbind(silhouette_average,silhouette_sd))
                                                             })
                                 })
    if(plot){
      pdf(paste0(name,'silhouette.pdf'),width=7,height=5)
      for(i in seq_along(K)){
        silhouette_average_plot=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) average_sd[order(average_sd[,1],decreasing=TRUE),1])),ncol=length(c))
        silhouette_sd_plot=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) average_sd[order(average_sd[,1],decreasing=TRUE),2])),ncol=length(c))
        silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
        matplot(matrix(rep(1:n_init,length(c))+rep(seq(-0.1,0.1,length.out=length(c)),each=n_init),ncol=length(c)),
                silhouette_average_plot,type='b',pch=16,lty=1,col=1+seq_along(c),ylim=c(0,1),xaxt='n',
                xlab='',ylab='Silhouette average',main=paste0('K=',K[i]))
        shift=seq(-0.1,0.1,length.out=length(c))
        for(ii in seq_along(c)){
          segments(1:n_init+shift[ii],silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
                   1:n_init+shift[ii],silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],col=ii+1)
          segments(1:n_init+shift[ii]-0.1,silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],
                   1:n_init+shift[ii]+0.1,silhouette_average_plot[,ii]+silhouette_sd_plot[,ii],col=ii+1)
          segments(1:n_init+shift[ii]-0.1,silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
                   1:n_init+shift[ii]+0.1,silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],col=ii+1)
          text(1:n_init+shift[ii],silhouette_average_plot[,ii]-silhouette_sd_plot[,ii],
               silhouette_order[,ii],col=ii+1,pos=1,cex=0.7)
        }
        legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),lty=1)
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
      }
      dev.off()
    }
    ### plot processing time ####################################################################################
    times=lapply(results,
                 function(results){
                   times=lapply(results,
                                function(results)
                                  times=unlist(lapply(results,function(results) results$time[3])))
                 })
    if(plot){
      pdf(paste0(name,'times.pdf'),width=7,height=5)
      y_max=max(unlist(times))*1.2
      times_plot=vector('list',length(K))
      for(i in seq_along(K)){
        times_plot[[i]]=Reduce(cbind,
                               lapply(seq_along(c),
                                      function(j)
                                        as.matrix(times[[i]][[j]][order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)]))
        )
        silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
        matplot(matrix(rep(1:n_init,length(c))+rep(seq(-0.1,0.1,length.out=length(c)),each=n_init),ncol=length(c)),
                times_plot[[i]],type='b',pch=16,lty=1,col=1+seq_along(c),ylim=c(0,y_max),xaxt='n',
                xlab='',ylab='Time',main=paste0('K=',K[i]))
        shift=seq(-0.1,0.1,length.out=length(c))
        for(ii in seq_along(c)){
          text(1:n_init+shift[ii],times_plot[[i]][,ii],silhouette_order[,ii],col=ii+1,pos=1,cex=0.7)
        }
        legend('topleft',legend=paste0('c=',c),col=1+seq_along(c),lty=1)
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
      }
      for(i in seq_along(K)){
        boxplot(times_plot[[i]],col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
                xlab='',ylab='Times',main=paste0('K=',K[i]))
      }
      dev.off()
    }
    
    ### plot dissimilarities ####################################################################################
    if(plot){
      D=lapply(results,
               function(results){
                 D=lapply(results,
                          function(results){
                            D=lapply(results,
                                     function(results){
                                       D=as.vector(results$probKMA_results$D)
                                     })
                          })
               })
      pdf(paste0(name,'dissimilarities.pdf'),width=7,height=5)
      y_max=max(unlist(D))
      for(i in seq_along(K)){
        D_plot=matrix(unlist(D[[i]]),ncol=length(c))
        boxplot(D_plot,col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
                xlab='',ylab='Dissimilarity',main=paste0('K=',K[i]))
      }
      for(i in seq_along(K)){
        for(j in seq_along(c)){
          silhouette_order=order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)
          D_plot=D[[i]][[j]][silhouette_order]
          boxplot(D_plot,col=j+1,ylim=c(0,y_max),names=silhouette_order,
                  xaxt='n',ylab='Dissimilarity',main=paste0('K=',K[i],'   c=',c[j]))
          axis(1,1:n_init,labels=rep('',n_init))
          mtext("Init",side=1,line=1)
        }
      }
      dev.off()
      D_clean=lapply(results,
                     function(results){
                       D_clean=lapply(results,
                                      function(results){
                                        D_clean=lapply(results,
                                                       function(results){
                                                         D_clean=as.vector(results$probKMA_results$D_clean)
                                                       })
                                      })
                     })
      pdf(paste0(name,'dissimilarities_clean.pdf'),width=7,height=5)
      y_max=max(unlist(D_clean))
      for(i in seq_along(K)){
        D_plot=matrix(unlist(D_clean[[i]]),ncol=length(c))
        boxplot(D_plot,col=1+seq_along(c),names=paste0('c=',c),ylim=c(0,y_max),
                xlab='',ylab='Dissimilarity',main=paste0('K=',K[i]))
      }
      for(i in seq_along(K)){
        for(j in seq_along(c)){
          silhouette_order=order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)
          D_plot=D_clean[[i]][[j]][order(silhouette_average_sd[[i]][[j]][,1],decreasing=TRUE)]
          boxplot(D_plot,col=j+1,ylim=c(0,y_max),names=silhouette_order,
                  xaxt='n',ylab='Dissimilarity',main=paste0('K=',K[i],'   c=',c[j]))
          axis(1,1:n_init,labels=rep('',n_init))
          mtext("Init",side=1,line=1)
        }
      }
      dev.off()
    }
    
    ### plot motif lengths ######################################################################################
    if(plot){
      motif_length=mapply(function(results){
        motif_length=mapply(function(results){
          motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                     function(results){
                                                       unlist(lapply(results$probKMA_results$V0,nrow))
                                                     })))
          return(as.matrix(motif_length))
        },results,SIMPLIFY=FALSE)
      },results,SIMPLIFY=FALSE)
      pdf(paste0(name,'lengths.pdf'),width=7,height=5)
      motif_length_plot=lapply(motif_length,
                               function(motif_length){
                                 lapply(motif_length,
                                        function(motif_length){
                                          motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                        })
                               })
      ymax=max(unlist(motif_length_plot))
      for(i in seq_along(K)){
        plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths',main=paste0('K=',K[i]))
        abline(h=c,col=1+seq_along(c))
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
        shift=seq(-0.1,0.1,length.out=length(c))
        silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
        for(ii in seq_along(c)){
          points(rep(1:n_init+shift[ii],each=K[i]),
                 motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
          text(rep(1:n_init+shift[ii],each=K[i]),motif_length[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
               rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
        }
        legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
      }
      dev.off()
      
      motif_clean_length=mapply(function(results){
        motif_length=mapply(function(results){
          motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                     function(results){
                                                       unlist(lapply(results$probKMA_results$V0_clean,nrow))
                                                     })))
          return(as.matrix(motif_length))
        },results,SIMPLIFY=FALSE)
      },results,SIMPLIFY=FALSE)
      pdf(paste0(name,'lengths_clean.pdf'),width=7,height=5)
      motif_length_plot=lapply(motif_clean_length,
                               function(motif_length){
                                 lapply(motif_length,
                                        function(motif_length){
                                          motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                        })
                               })
      ymax=max(unlist(motif_length_plot))
      for(i in seq_along(K)){
        plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths',main=paste0('K=',K[i]))
        abline(h=c,col=1+seq_along(c))
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
        shift=seq(-0.1,0.1,length.out=length(c))
        silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
        for(ii in seq_along(c)){
          points(rep(1:n_init+shift[ii],each=K[i]),
                 motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
          text(rep(1:n_init+shift[ii],each=K[i]),motif_clean_length[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
               rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
        }
        legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
      }
      dev.off()
      
      pdf(paste0(name,'lengths_perc.pdf'),width=7,height=5)
      motif_length_perc=mapply(function(results){
        motif_length=mapply(function(results){
          motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                     function(results){
                                                       unlist(lapply(results$probKMA_results$V0,nrow))/results$probKMA_results$c*100
                                                     })))
          return(as.matrix(motif_length))
        },results,SIMPLIFY=FALSE)
      },results,SIMPLIFY=FALSE)
      motif_length_plot=lapply(motif_length_perc,
                               function(motif_length){
                                 lapply(motif_length,
                                        function(motif_length){
                                          motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=0.05),ncol=n_init)
                                        })
                               })
      ymax=max(unlist(motif_length_plot))
      for(i in seq_along(K)){
        plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths (% of minimum length)',main=paste0('K=',K[i]))
        abline(h=100)
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
        shift=seq(-0.1,0.1,length.out=length(c))
        silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
        for(ii in seq_along(c)){
          points(rep(1:n_init+shift[ii],each=K[i]),
                 motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
          text(rep(1:n_init+shift[ii],each=K[i]),motif_length_perc[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
               rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
        }
        legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
      }
      dev.off()
      
      pdf(paste0(name,'lengths_clean_perc.pdf'),width=7,height=5)
      motif_clean_length_perc=mapply(function(results){
        motif_length=mapply(function(results){
          motif_length=as.matrix(Reduce(cbind,lapply(results,
                                                     function(results){
                                                       unlist(lapply(results$probKMA_results$V0_clean,nrow))/results$probKMA_results$c*100
                                                     })))
          return(as.matrix(motif_length))
        },results,SIMPLIFY=FALSE)
      },results,SIMPLIFY=FALSE)
      motif_length_plot=lapply(motif_clean_length_perc,
                               function(motif_length){
                                 lapply(motif_length,
                                        function(motif_length){
                                          motif_length=motif_length+matrix(rnorm(nrow(motif_length)*n_init,mean=0,sd=10^-1),ncol=n_init)
                                        })
                               })
      ymax=max(unlist(motif_length_plot))
      for(i in seq_along(K)){
        plot(1:n_init,type="n",xaxt='n',xlab='',ylim=c(-1,ymax),ylab='Motif lengths (% of minimum length)',main=paste0('K=',K[i]))
        abline(h=100)
        axis(1,1:n_init,labels=rep('',n_init))
        mtext("Init",side=1,line=1)
        shift=seq(-0.1,0.1,length.out=length(c))
        silhouette_order=matrix(unlist(lapply(silhouette_average_sd[[i]],function(average_sd) order(average_sd[,1],decreasing=TRUE))),ncol=length(c))
        for(ii in seq_along(c)){
          points(rep(1:n_init+shift[ii],each=K[i]),
                 motif_length_plot[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],col=ii+1,pch=8)
          text(rep(1:n_init+shift[ii],each=K[i]),motif_clean_length_perc[[i]][[ii]][,order(silhouette_average_sd[[i]][[ii]][,1],decreasing=TRUE)],
               rep(silhouette_order[,ii],each=K[i]),col=ii+1,pos=1,cex=0.7)
        }
        legend('bottomleft',legend=paste0('c=',c),col=1+seq_along(c),pch=8)
      }
      dev.off()
    }
    
    ### output ##################################################################################################
    find_candidate_motifs_results <- list(name=name,K=K,c=c,n_init=n_init,silhouette_average_sd=silhouette_average_sd,times=times)
    
    save(find_candidate_motifs_results, file = paste0(name,'find_candidate_motifs_results.RData'))
    }
    ### filter candidate motifs based on silhouette average and size
    silhouette_average = Reduce(rbind, Reduce(rbind, find_candidate_motifs_results$silhouette_average_sd))[ , 1] # retrieve silhouette average for all candidate motifs
    pb$update(0.5)
    filter_candidate_motifs_results = filter_candidate_motifs(find_candidate_motifs_results,
                                                              sil_threshold = quantile(silhouette_average, probKMA_options$sil_threshold),
                                                              size_threshold = 2)
    ### cluster candidate motifs based on their distance and select radii
    pb$update(0.75)
    cluster_candidate_motifs_results = cluster_candidate_motifs(filter_candidate_motifs_results,
                                                                motif_overlap = 0.6,
                                                                worker_number = worker_number)
    cluster_candidate_motifs_results$transformed = probKMA_options$transformed
    ### plot cluster candidate motifs results
    pdf(paste0(name,'FMD_clustering_candidate_motifs.pdf'), height = 12, width = 9)
    cluster_candidate_motifs_plot(cluster_candidate_motifs_results, ask = FALSE)
    dev.off()
    ### search selected motifs
    cluster_candidate_motifs_results$Y0 <- Y0
    cluster_candidate_motifs_results$Y1 <- Y1
    
    pb$update(0.9)
    motifs_search_results = motifs_search(cluster_candidate_motifs_results,
                                          use_real_occurrences = FALSE, length_diff = +Inf,
                                          worker_number = worker_number)
    
    ### plot FMD results (NB: no threshold of frequencies of motif found!)
    pdf(paste0(name,'FMD_results.pdf'), height = 7, width = 17)
    motifs_search_plot(motifs_search_results, ylab = 'x(t)', freq_threshold = 1,
                       transformed = probKMA_options$transformed)
    dev.off()
    
    save(find_candidate_motifs_results, silhouette_average, filter_candidate_motifs_results,
         cluster_candidate_motifs_results, motifs_search_results,
         file=paste0(name,'FMD_results.RData'))
    pb$update(1.0)
    return(list(find_candidate_motifs_results = find_candidate_motifs_results,
                silhouette_average = silhouette_average,
                filter_candidate_motifs_results = filter_candidate_motifs_results,
                cluster_candidate_motifs_results = cluster_candidate_motifs_results,
                motifs_search_results = motifs_search_results))
  }else if(method=="FunBIalign")
  {
    if(length(funBIalign_options) != 3) {
      stop('\'funBIalign_options\' must contain: \'portion_len\',\'min_card\',\'cut_off\'.')
    }
    portion_len <- funBIalign_options$portion_len
    min_card <- funBIalign_options$min_card
    cut_off <- funBIalign_options$cut_off
    # check required input
    if(missing(portion_len) || is.null(portion_len))
      stop('portion_len attribute must be specified')
    if(missing(min_card) || is.null(min_card))
      stop('min_card attribute must be specified')
    
    pb <- progress_bar$new(
      format = "Progress [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = 5,  
      width = 50
    )
    
    pb$update(0.1)
    #compute maximum length
    maxLen <- max(sapply(Y0, nrow))
    full_data <- sapply(Y0, padding, maxLen,simplify = FALSE)
    
    # transform list into a matrix where each curve is on a row
    full_data <- t(sapply(full_data,cbind))

    window_data <- NULL
    list_of_recommendations_ordered <- NULL
    vec_of_scores_ordered <- NULL
    if(!file.exists(paste0(name,'/resFunBi.rds')))
    {
      dir.create(paste0(name),showWarnings=TRUE)
  cppFunction('
  Rcpp::List createWindow(Rcpp::NumericMatrix& data,
                          unsigned int portion_len,
                          unsigned int worker_number){
  
  const unsigned int totdim = data.ncol();
  const unsigned int totobs = data.nrow();
  const unsigned int totrows = (totdim - portion_len + 1) * totobs;
  const unsigned int totportion = totdim - portion_len + 1;
  
  // set the size for data structure
  const arma::mat dataRef(data.begin(), totobs, totdim,false,true);
  Rcpp::NumericMatrix windowData(totrows,portion_len);
  arma::mat windowDataRef(windowData.begin(),totrows,portion_len,false,true);
  std::vector<std::string> window_rownames(totrows);
  
  if(totobs <= 1) // we have a single curve
  {
    // fill in my data 
#ifdef _OPENMP
#pragma omp parallel for num_threads(worker_number)
#endif
    for(arma::uword i = 0; i < totrows; ++i)
    {
      windowDataRef.row(i) = dataRef.cols(i,i + portion_len - 1); //Each peace of curve is stored in a row
      window_rownames[i] = "1_" + std::to_string(i + 1) + "_" + std::to_string(i + portion_len);         
    }
  }
  else // we have multiple curves
  {
    // fill in my data
#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(worker_number)
#endif
      for(arma::uword k = 0; k < totobs; ++k)
      {
        for (arma::uword i = 0; i < totportion; ++i)
        {
          windowDataRef.row(totportion * k + i) = dataRef(k,arma::span(i,i + portion_len - 1));
          window_rownames[i + k * totportion] = std::to_string(k+1) + "_" + std::to_string(i + 1) + "_" + std::to_string(i + portion_len); 
        }
                
      }
  }
  return Rcpp::List::create(windowData,window_rownames);
}',depends = "RcppArmadillo")
    # step 1
    window_data_list <- createWindow(full_data,portion_len,worker_number)
    window_data <- window_data_list[[1]]
    rownames(window_data) <- window_data_list[[2]]
    pb$update(0.3)
    # set it as a matrix and omit rows with NAs
    window_data <- na.omit(window_data)
    
    # compute numerosity 
    numerosity <- rep(dim(full_data)[2] - portion_len + 1,times = dim(full_data)[1])
    removed_rows <- setdiff(window_data_list[[2]],rownames(window_data))
    curve_indices <- as.numeric(gsub("^(\\d+)_.*", "\\1", removed_rows))
    tab <- table(curve_indices)
    numerosity[as.numeric(names(tab))] <- numerosity[as.numeric(names(tab))] - tab
    
    rm(window_data_list)
    rm(removed_rows)
    rm(tab)
    
  cppFunction('
  Rcpp::NumericMatrix createDistance(Rcpp::NumericMatrix& windowData,
                                     const arma::vec& numerosity,
                                     unsigned int worker_number){
  
  int outrows = windowData.nrow();
  int outcols = windowData.ncol();
  bool isSingle = (numerosity.size() == 1);
  
  double outcols_inv = 1.0 / static_cast<double>(outcols);
  const arma::mat windowDataRef(windowData.begin(),outrows,outcols,false,true);
  const arma::colvec& vsum = arma::sum(windowDataRef,1) * outcols_inv;
  Rcpp::NumericMatrix result(outrows,outrows);
  arma::mat scoreData(result.begin(),outrows,outrows,false,true);
  
#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(worker_number)
#endif
  for (arma::uword i = 1; i < outrows; ++i) { // for any functional observation
    const arma::rowvec& x_i = windowDataRef.row(i); // curve i 
    for (arma::uword j = 0; j < i; ++j) {
      
      // Precompute common terms
      const arma::rowvec& crossSum = x_i + windowDataRef.row(j);;
      const double commonTerm = arma::accu(crossSum) * (outcols_inv * 0.5);
      
      // fill the lower part
      scoreData(i,j) = arma::accu(arma::square(
                                   x_i - vsum[i] - crossSum/2
                                   + commonTerm))*outcols_inv;
    }
  }
  
  // Modify accolites 
  double M = static_cast<double>(*std::max_element(scoreData.begin(),scoreData.end())) + 1000.0; // very large distance for accolites
  int overlap = static_cast<int>(std::floor(outcols / 2.0)); // number of right/left accolites
  
  // Single curve
  if (isSingle) {
    for (int i = 0; i < outrows-1; ++i) {
      // check for right accolites
      arma::uword start = i+1;
      arma::uword end = std::min<int>(i+overlap,outrows-1);
      scoreData(arma::span(start,end),i) += M;
    }
    
    // Multiple curves
  } else {
    int until_here = 0;
    int numerosity_size = numerosity.size();
    for(int j = 0; j <numerosity_size;++j)
    {
      int num = numerosity[j] - 1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(worker_number)
#endif
      for(int i = until_here; i < until_here + numerosity[j]-1;++i)
      {
        arma::uword start = i+1;
        arma::uword end = std::min<int>(std::min<int>(i+overlap ,i + (num--)),outrows-1);
        scoreData(arma::span(start,end),i) += M;
      }
      until_here += numerosity[j];
    }
  }
  return Rcpp::wrap(result);
}',depends = "RcppArmadillo")
# step 2: compute fMRS-based dissimilarity matrix
    D_fmsr <- createDistance(window_data,numerosity,worker_number)
    pb$update(0.5)
    rownames(D_fmsr) <- colnames(D_fmsr) <- rownames(window_data)
    ## STEP 3 -----
    # step 3: get the sub-trees (tree_s)
    minidend  <- get_minidend(as.dist(D_fmsr))
    # step 3: identify seeds and corresponding families

    all_paths   <- get_path_complete(minidend, window_data, min_card = min_card,worker_number = worker_number)
    all_paths   <- all_paths[!(lapply(all_paths, is.null) %>% unlist())] 
    pb$update(0.65)
    # step 3: get recommended nodes and their info (cardinality and score)
    # collect all recommended nodes (as an array of portion ids)
    all_recommended_labels <- lapply(all_paths, function(x){x$recommended_node_labels})
    list_of_recommendations <- lapply(rapply(
      all_recommended_labels, enquote, how='unlist'),
      eval)
    # get recommended node cardinality
    vec_of_card <- lapply(list_of_recommendations, length) %>% unlist()
    
    # get vector of adjusted fMSR
    vec_of_scores <- lapply(all_paths, function(x){x$recommended_node_scores}) %>% unlist()
    
    resFunBiMeta = list(window_data = window_data,
                        list_of_recommendations = list_of_recommendations,
                        vec_of_scores = vec_of_scores)
    # save results 
    saveRDS(resFunBiMeta,file = paste0(name,'/resFunBi.rds')) 
    
    } else {
      resFunBiMeta <- readRDS(paste0(name,'/resFunBi.rds'))
      window_data <- resFunBiMeta$window_data
      list_of_recommendations <-  resFunBiMeta$list_of_recommendations
      vec_of_scores <- resFunBiMeta$vec_of_scores
    }
    
    ## STEP 4 -----
    # STEP 4: post-processing and rearranging results using different criteria
    pb$update(0.8)
    if (stopCriterion == "Variance") {
      motif_var <- lapply(list_of_recommendations, 
                          function(x){
                            temp <- window_data[x,]
                            temp_mean <- temp %>% colMeans()
                            ((temp - temp_mean)^2 %>% sum())/(nrow(temp))
                          }) %>% unlist()
      
      var_order <- motif_var %>% order(decreasing = TRUE)
      vec_of_scores_ordered <- vec_of_scores[var_order]
      list_of_recommendations_ordered <- list_of_recommendations[var_order]
    } else {
      ### CRITERION: adjusted fMSR ----
      best_order <- vec_of_scores %>% order()
      vec_of_scores_ordered <- vec_of_scores[best_order]
      # ordered list_of_recommendations
      list_of_recommendations_ordered <- list_of_recommendations[best_order]  
    }
    
    #for every recommended motif, compute all its accolites
    all_accolites <- lapply(list_of_recommendations_ordered,
                            function(x){
                              lapply(x, get_accolites, window_data, portion_len, FALSE) %>% unlist()
                            })
    
    # Starting from the top, we compare each motif to those with higher rank. If all portions of
    # are acolytes to portions of an higher ranking motif, we filter it out; 
    # otherwise we retain it.

    # we identify the ones to delete
    cppFunction('Rcpp::IntegerVector deleteV(const Rcpp::List& list_of_recommendations_ordered,
                                             const Rcpp::List&  all_accolites,
                                             const Rcpp::NumericVector& vec_of_scores_ordered){
  int n = list_of_recommendations_ordered.size();
  Rcpp::IntegerVector del;
  
  for (int i = 1; i < n; ++i) { // compare every motif (from the second one)
    const Rcpp::CharacterVector& node_1 = list_of_recommendations_ordered[i];
    for (int j = 0; j < i; ++j) { // to the other higher ranked nodes
      const Rcpp::CharacterVector& node_2 = list_of_recommendations_ordered[j];
      if (node_1.size() <= node_2.size()) {
        const Rcpp::CharacterVector& accolites_2 = all_accolites[j];
        if (vec_of_scores_ordered(i) > vec_of_scores_ordered(j) && 
            Rcpp::is_true(Rcpp::all(Rcpp::in(node_1, accolites_2)))){
          del.push_back(i+1);
          break;
        }
      }
    }
  }
  
  return del;
  
}',depends="RcppArmadillo")
    delete <- deleteV(list_of_recommendations_ordered,
                      all_accolites,
                      vec_of_scores_ordered)
    # we delete the recommended nodes and we order the remaining ones
    list_of_recommendations_ordered <- list_of_recommendations_ordered[-delete]
    vec_of_scores_ordered <- vec_of_scores_ordered[-delete]
    if(!is.null(cut_off)) {
      valid_pos <- which(vec_of_scores_ordered < cut_off)
      vec_of_scores_ordered <- vec_of_scores_ordered[valid_pos]
      list_of_recommendations_ordered <- list_of_recommendations_ordered[valid_pos] 
      if(is.null(vec_of_scores_ordered))
        stop("the value of \'cut_off\' is to high. No motifs discovered. Provide a lower value")
    }
    pb$update(0.9)
    # Creating some plot ----
    ## Plot the data and highlight the motif occurrences in red -----
    if(plot)
    {
      pdf(paste0(name,"/plot_",stopCriterion,".pdf"),width=10,height=5)
      layout(matrix(1:2,ncol=2,byrow=TRUE),widths=c(8.5,1))
      for(q in 1:length(list_of_recommendations_ordered)){
        temp_motif <- list_of_recommendations_ordered[[q]]
        lots_in_motif <- lapply(temp_motif,
                                function(x){
                                  strsplit(x, '_') %>% 
                                    unlist() %>%
                                    as.numeric()}) %>% 
          unlist() %>%
          matrix(ncol=3, byrow=T)
        current_curves <- unique(as.numeric(gsub("[^0-9]", "", substr(temp_motif, 1, 2))))
        title   <- paste0("Number of instances: ", length(temp_motif),
                          " - adj fMSR:  ", vec_of_scores_ordered[q] %>% round(3))
        par(mar = c(5, 4, 2, 2) + 0.1)
        matplot(t(full_data), ylab='', xlab='',
                lwd=1.5,lty = 1, type = 'l', main = title,col=alpha('gray30',0.15))
        matplot(matrix(full_data[current_curves,],ncol=length(current_curves),byrow = TRUE), ylab='', xlab='',
                lwd=1.5,lty = 5, type = 'l', main = title,col=rainbow(length(current_curves)),add=TRUE)
        rect(lots_in_motif[,2], min(full_data,na.rm = TRUE)-10, lots_in_motif[,3], max(full_data,na.rm = TRUE) +10,
             border = alpha("firebrick3", 0.05), col = alpha("firebrick3", 0.05))
      
        lapply(1:nrow(lots_in_motif), function(k) {
          matplot(lots_in_motif[k, 2]:lots_in_motif[k, 3],
                  full_data[lots_in_motif[k, 1], lots_in_motif[k, 2]:lots_in_motif[k, 3]],
                  type = 'l', add = TRUE, col = 'red', lwd = 2.5)
          #box(col = "grey40", lwd = 2)
          #axis(1, at = seq(1, length(Y0), by = 12), col = "grey40", col.ticks = "grey40", col.axis = "grey60", cex.axis = 1.5)
          #axis(2, col = "grey40", col.ticks = "grey40", col.axis = "grey60", cex.axis = 1.5)
        })
        par(mar=c(0,0,0,0))
        plot.new()
        legend("topright", 
               legend = c(paste0("c",current_curves),paste0('motif_',q)),
               col = c(rainbow(length(current_curves)),'red'),lwd=c(rep(1,length(current_curves)),4), lty = 1,cex=0.75)
      }
      layout(matrix(1:2,ncol=2,byrow=TRUE),widths=c(5,1))
      ## Plot only the motif -----
      for(q in 1:length(list_of_recommendations_ordered)){
        par(mar = c(5, 4, 2, 2) + 0.1)
        temp_motif <- list_of_recommendations_ordered[[q]]
        current_curves <- as.numeric(gsub("[^0-9]", "", substr(temp_motif, 1, 2)))
        plot_me <- t(window_data[temp_motif,])
        motifMean <- rowMeans(plot_me)
        title   <- paste0("Number of occurrences: ", dim(plot_me)[2],
                          " - adj fMSR:  ", vec_of_scores_ordered[q] %>% round(3))
        
        matplot(plot_me, ylab='', xlab='',type='l',
                col=seq_len(length(temp_motif))+1,lwd=2,lty=5, main = title)
        lines(motifMean,col='black',lwd=5,lty=1)
        par(mar=c(0,0,0,0))
        plot.new()
        legend('topright',
               legend=c(paste0('c',current_curves),'motif center'),
               col=c(seq_len(length(temp_motif))+1,'black'),lwd=c(rep(2,length(temp_motif)),5),lty=c(rep(5,length(temp_motif)),1),
               xpd=TRUE,cex=0.75)
        #box(col="grey40", lwd = 2)
        #axis(1, at=1:portion_len, labels = 1:portion_len, col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
        #axis(2, col="grey40", col.ticks="grey40", col.axis="grey60", cex.axis=1.5)
      }
      dev.off()
    }
    pb$update(1.0)
    resFunBi = list(list_of_recommendations_ordered = list_of_recommendations_ordered,
                    vec_of_scores_ordered = vec_of_scores_ordered)
    return(resFunBi)
  } else {
    stop('\'method\' not found. It must be choosen between \'probKMA\' and \'funBIalign\'')
  }
}