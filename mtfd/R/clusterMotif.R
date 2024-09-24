#' @title clusterMotif
#'
#' @description Run multiple times probKMA function with different K,c and initializations,
#' with the aim to find a set of candidate motifs.
#' If the folder name_KK_cc is already present and n result files are already present,
#' load them and continue with the n_init-n runs.
#'
#' @param Y0 list of N vectors, for univariate curves y_i(x), or list of N matrices with d columns,
#' for d-dimensional curves y_i(x),  with the evaluation of curves (all curves should be evaluated
#' on a uniform grid). When y_j(x)=NA in the dimension j, then y_j(x)=NA in ALL dimensions
#' @param Y1 list of N vectors, for univariate derivative curves y'_i(x), or
#' list of N matrices with d columns, for d-dimensional derivative curves y'_i(x),
#' with the evaluation of the curves derivatives (all curves should be evaluated on a uniform grid).
#' When y'_j(x)=NA in the dimension j, then y'_j(x)=NA in ALL dimensions.
#' Must be provided when diss='d1_L2' or diss='d0_d1_L2'.
#' @param K vector with numbers of motifs that must be tested.
#' @param c vector with minimum motifs lengths that must be tested.
#' @param n_init number of random initialization for each combination of K and c.
#' @param name name of the folders when the results are saved.
#' @param names_var vector of length d, with names of the variables in the different dimensions.
#' @param probKMA_options list with options for probKMA (see the help of probKMA).
#' @param silhouette_align True or False. If True, try all possible alignments between the curve pieces
#' @param criterium reorder the results based on fMRS or Variance.
#' when calculating the adapted silhouette index on the results of probKMA
#' @param plot if TRUE, summary plots are drawn.
#' @return A list containing: K, c, n_init and name;...
#' @return \item{times}{ list of execution times of ProbKMA for each combination of K, c, and n_init}
#' @return \item{silhouette_average_sd}{ list of the mean (silhouette_average) and standard deviation (silhouette_sd) of the silhouette indices for each execution of the ProbKMA function}
#' @author Niccolò Feresini, Riccardo Lazzarini
#' @export

clusterMotif <- function(Y0,method,portion_len = NULL,min_card = NULL,criterium="fMRS",cut_off=NULL,
                         Y1=NULL,P0=matrix(),S0=matrix(),plot=TRUE,K=NULL,c=NULL,n_init=10,diss="d0_L2",alpha=0,
                         stopCriterion="max",return_options=TRUE,name='results',names_var='',
                         probKMA_options=list(standardize=TRUE,c_max=Inf,
                                              iter_max=1e3,iter4elong=10,
                                              trials_elong=10,max_gap=0.2,
                                              quantile=0.25,tol=1e-8,
                                              tol4elong=1e-3,max_elong=0.5,
                                              deltaJK_elong=0.05,iter4clean=50,
                                              tol4clean=1e-4,quantile4clean=0.5,
                                              m=2,w=1,transformed = FALSE,silhouette_align=FALSE,
                                              n_subcurves = 10,sil_threshold=0.9),
                         seed = 1,exe_print = FALSE,set_seed = FALSE,
                         worker_number = NULL,
                         V_init=NULL, n_init_motif = 0){
  
  
  if(method == 'ProbKMA')
  {
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
    }
    if(is.null(worker_number))
      worker_number <- core_number-1
    
    rm(core_number)
    len_mean=mean(unlist(lapply(Y0,nrow)))
    c_mean=mean(c)
    
    arguments = list(Y0 = Y0, 
                     Y1 = Y1,
                     P0 = P0, 
                     S0 = S0,
                     standardize=probKMA_options$standardize,
                     c_max = probKMA_options$c_max,
                     iter_max = probKMA_options$iter_max,
                     iter4elong = probKMA_options$iter4elong,
                     max_gap = probKMA_options$max_gap,
                     quantile = probKMA_options$quantile, 
                     tol = probKMA_options$tol,
                     trials_elong = probKMA_options$trials_elong,
                     tol4elong = probKMA_options$tol4elong, 
                     max_elong = probKMA_options$max_elong, 
                     deltaJK_elong = probKMA_options$deltaJK_elong, 
                     iter4clean = probKMA_options$iter4clean, 
                     tol4clean = probKMA_options$tol4clean,
                     quantile4clean = probKMA_options$quantile4clean,
                     m = probKMA_options$m,
                     w = probKMA_options$w, 
                     transformed = probKMA_options$transformed,
                     return_options = return_options,
                     alpha = alpha,
                     diss = diss,
                     stopCriterion = stopCriterion, 
                     seed = seed, 
                     K = NULL, #any value
                     c = NULL, #any value
                     v_init = NULL, #any value
                     exe_print = exe_print,
                     set_seed = set_seed,
                     align = probKMA_options$silhouette_align,
                     sil_threshold = probKMA_options$sil_threshold,
                     n_threads = worker_number)
    if(((len_mean/c_mean>10)&(length(K)*length(c)*n_init<40))|(length(K)*length(c)*n_init==1)){
      if(is.null(probKMA_options$worker_number))
        probKMA_options$worker_number=worker_number
      cl_find=NULL
    }else{
      probKMA_options$worker_number=1
      if(worker_number>1){
        cl_find=parallel::makeCluster(worker_number,timeout=60*60*24*30)
        parallel::clusterExport(cl_find,c('name','names_var',
                                          'probKMA_plot','probKMA_silhouette_plot',
                                          '.mapply_custom','.diss_d0_d1_L2','.domain',
                                          '.select_domain','.find_min_diss',
                                          'probKMA_wrap','arguments'),envir=environment()) 
        parallel::clusterCall(cl_find,function()library(parallel,combinat))
        on.exit(parallel::stopCluster(cl_find))
      }else{
        cl_find=NULL
      }
    }
    
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
            results = do.call(mtfd:::probKMA_wrap,arguments)
            probKMA_results = results[[1]]
            silhouette_results = results[[2]]
            end=proc.time()
            time=end-start
            iter=probKMA_results$iter
            iter_max=arguments$iter_max
            if(iter==iter_max)
              warning('Maximum number of iteration reached. Re-starting.')
          }
          transformed=probKMA_results$transformed
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
          mtfd:::probKMA_plot(probKMA_results,ylab=names_var,cleaned=FALSE,transformed=transformed) #transformed
          dev.off()
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
          mtfd:::probKMA_plot(probKMA_results,ylab=names_var,sil_avg = silhouette_results[[4]],cleaned=TRUE,transformed=transformed)
          dev.off()
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
          silhouette=mtfd:::probKMA_silhouette_plot(results, 
                                             align=arguments$silhouette_align,
                                             plot=TRUE)
          dev.off()
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
          if(N==1) {
            cat("Number of curves is equal to \'1\'. Applying random split.")
            silh_value = -1.0
            # Function to compute variance in a sliding window, ignoring completely empty/NA curves
            .compute_local_variance <- function(curve, window_size) {
              apply(curve, 2, function(dim_curve) {
                zoo::rollapply(dim_curve, width = window_size, FUN = function(x) {
                  if (all(is.na(x))) {
                    return(NA)
                  } else {
                    return(var(na.omit(x)))
                  }
                }, fill = NA, align = "center", partial = TRUE)
              })
            }
            
            # Function to compute the slope (first derivative) in a sliding window for each dimension
            .compute_local_slope <- function(curve, window_size) {
              apply(curve, 2, function(dim_curve) {
                zoo::rollapply(dim_curve, width = window_size, FUN = function(x) {
                  if (length(na.omit(x)) > 1) {
                    return(lm(na.omit(x) ~ seq_along(na.omit(x)))$coefficients[2])
                  } else {
                    return(NA)
                  }
                }, fill = NA, align = "center", partial = TRUE)
              })
            }
            
            .split_curve_statistically <- function(curve, expected_motif_length, K, min_length, window_size, n_subcurves,randomness_factor = 0.3) {
              length_curve <- nrow(curve)  # Number of observations in the curve
              
              # Step 1: Compute local variance and slope to detect changes in the curve
              local_variance <- .compute_local_variance(curve, window_size)
              local_slope <- .compute_local_slope(curve, window_size)
              
              # Normalize the variance and slope for fair comparison
              normalized_variance <- scale(local_variance, center = TRUE, scale = TRUE)
              normalized_slope <- scale(local_slope, center = TRUE, scale = TRUE)
              
              # Combine variance and slope into a single score
              combined_score <- apply(normalized_variance + normalized_slope, 1, function(x) mean(x, na.rm = TRUE))
              
              # Step 2: Identify candidate split points where the score is high (indicating a transition)
              threshold <- quantile(combined_score, 0.75, na.rm = TRUE)  # Top 25% of the score
              candidate_splits <- which(combined_score >= threshold & !is.na(combined_score))
              
              # Step 3: Create valid split points while ensuring each valid portion of values meets the expected length requirement
              split_points <- c(1)  # Start point
              last_split <- 1
              for (i in candidate_splits) {
                # Check if we can add a split while respecting minimum length
                if (i - last_split >= min_length) {
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
              
              # Add last segment if it meets the minimum length requirement
              if (length_curve - last_split >= min_length) {
                split_points <- c(split_points, length_curve)
              } else
              {
                split_points[length(split_points)] <- length_curve
              }
              
              split_points <- unique(split_points)
              
              # Step 4: Adjust the number of split points to match n_subcurves if provided
              if (length(split_points) - 1 > n_subcurves) {
                # Too many splits, reduce by merging closest split points
                while (length(split_points) - 1 > n_subcurves) {
                  # Find the smallest gap between consecutive split points
                  gaps <- diff(split_points)
                  min_gap_index <- which.min(gaps)
                  split_points <- split_points[-(min_gap_index + 1)]  # Remove the smaller split
                }
              } else if (length(split_points) - 1 < n_subcurves) {
                # Too few splits, add new split points in the largest gaps between existing points
                additional_splits_needed <- n_subcurves - (length(split_points) - 1)
                local_iter <- 0
                while (additional_splits_needed > 0) {
                  local_iter <- local_iter + 1
                  # Find the largest gap between consecutive split points
                  gaps <- diff(split_points)
                  max_gap_index <- which.max(gaps)
                  # Add a split point in the middle of the largest gap, respecting min_length
                  new_split <- round(mean(split_points[max_gap_index:(max_gap_index + 1)]))
                  
                  # Ensure the new split doesn't violate minimum length
                  if (new_split - split_points[max_gap_index] >= min_length && 
                      split_points[max_gap_index + 1] - new_split >= min_length) {
                    split_points <- sort(c(split_points, new_split))
                    additional_splits_needed <- additional_splits_needed - 1
                  } else {
                    stop(paste0("Failed to find valid split points. Unable to divide the curve into subcurves where each subcurve has a minimum length of ", c, " observations. Please check the curve data or adjust the minimum length constraint."))
                  }
                }
              }
              
              # Step 5: Create subcurves based on valid split points
              subcurves <- list()
              for (i in 1:(length(split_points) - 1)) {
                subcurve <- curve[split_points[i]:split_points[i + 1], ]
                if (nrow(as.matrix(subcurve)) > 0) {  # Ensure subcurve is not empty
                  subcurves[[i]] <- subcurve
                }
              }
              return(list(subcurves = subcurves, split_points = split_points))
            }
            counter <- 0
            while(iter==iter_max && !all(silh_value > arguments$sil_threshold)){
              counter <- counter + 1
              start=proc.time()
              small_seed = small_seed + 1
              if(arguments$diss == 'd0_d1_L2') {
                split_results <- .split_curve_statistically(Y0[[1]], c, K, c, c,probKMA_options$n_subcurves, randomness_factor = 0.0)
                arguments$Y0 <- split_results$subcurves
                Y1_subcurves <- list()
                for (i in seq_along(split_results$split_points)[-length(split_results$split_points)]) {
                  Y1_subcurves[[i]] <- Y1[[1]][pos[i]:pos[i + 1]]
                }
                arguments$Y1 <- Y1_subcurves
              }
              else if(arguments$diss == 'd0_L2') {
                split_results <- .split_curve_statistically(Y0[[1]], c, K, c, c,probKMA_options$n_subcurves, randomness_factor = 0.0)
                arguments$Y0 <- split_results$subcurves
              }
              else if(arguments$diss == 'd1_L2') {
                split_results <- .split_curve_statistically(Y1[[1]], c, K, c, c,probKMA_options$n_subcurves,randomness_factor = 0.0)
                arguments$Y1 <- split_results$subcurves
              }
              arguments$K = K
              arguments$c = c
              arguments$quantile4clean = 1/K
              results = do.call(mtfd:::probKMA_wrap,arguments)
              probKMA_results = results[[1]]
              silhouette_results = results[[2]]
              silh_value = silhouette_results$silhouette_average
              end=proc.time()
              time=end-start
              iter=probKMA_results$iter
              iter_max=arguments$iter_max
              
              if(iter==iter_max)
                warning('Maximum number of iteration reached. Re-starting.')
              
              if(counter == 1000)
                stop('No available splits found in 1000 trials. Please decrease the \'sil_threshold\'.')
            }
          } else {
            while(iter==iter_max){
              start=proc.time()
              small_seed = small_seed + 1
              arguments$K = K
              arguments$c = c
              arguments$quantile4clean = 1/K
              results = do.call(mtfd:::probKMA_wrap,arguments)
              probKMA_results = results[[1]]
              silhouette_results = results[[2]]
              end=proc.time()
              time=end-start
              iter=probKMA_results$iter
              iter_max=arguments$iter_max
              
              if(iter==iter_max)
                warning('Maximum number of iteration reached. Re-starting.')
            }
          }
          
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'.pdf'),width=20,height=10)
          mtfd:::probKMA_plot(probKMA_results,
                              sil_avg = silhouette_results[[4]],
                              ylab=names_var,
                              cleaned=FALSE) 
          dev.off()
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'clean.pdf'),width=20,height=10)
          mtfd:::probKMA_plot(probKMA_results,
                              sil_avg = silhouette_results[[4]],
                              ylab=names_var,
                              cleaned=TRUE) 
          dev.off()
      
          pdf(paste0(name,"K",K,"_c",c,'/random',i,'silhouette.pdf'),width=7,height=10)
          silhouette = mtfd:::probKMA_silhouette_plot(silhouette_results,K,plot = plot)
          dev.off()
          
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
    filter_candidate_motifs_results = mtfd:::filter_candidate_motifs(find_candidate_motifs_results,
                                                              sil_threshold = quantile(silhouette_average, probKMA_options$sil_threshold),
                                                              size_threshold = 2)
    
    ### cluster candidate motifs based on their distance and select radii
    cluster_candidate_motifs_results = mtfd:::cluster_candidate_motifs(filter_candidate_motifs_results,
                                                                       motif_overlap = 0.6)
    cluster_candidate_motifs_results$transformed = probKMA_options$transformed
    ### plot cluster candidate motifs results
    pdf(paste0(name,'FMD_clustering_candidate_motifs.pdf'), height = 12, width = 9)
    mtfd:::cluster_candidate_motifs_plot(cluster_candidate_motifs_results, ask = FALSE)
    dev.off()
    ### search selected motifs
    cluster_candidate_motifs_results$Y0 <- Y0
    cluster_candidate_motifs_results$Y1 <- Y1
    motifs_search_results = mtfd:::motifs_search(cluster_candidate_motifs_results,
                                          use_real_occurrences = FALSE, length_diff = +Inf)
    
    ### plot FMD results (NB: no threshold of frequencies of motif found!)
    pdf(paste0(name,'FMD_results.pdf'), height = 7, width = 17)
    mtfd:::motifs_search_plot(motifs_search_results, ylab = 'x(t)', freq_threshold = 1)
    dev.off()
    
    save(find_candidate_motifs_results, silhouette_average, filter_candidate_motifs_results,
         cluster_candidate_motifs_results, motifs_search_results,
         file=paste0(name,'FMD_results.RData'))
    return(list(find_candidate_motifs_results = find_candidate_motifs_results,
                silhouette_average = silhouette_average,
                filter_candidate_motifs_results = filter_candidate_motifs_results,
                cluster_candidate_motifs_results = cluster_candidate_motifs_results,
                motifs_search_results = motifs_search_results))
  }
  if(method=="FunBialign")
  {
    # check required input
    if(missing(portion_len) || is.null(portion_len))
      stop('portion_len attribute must be specified')
    if(missing(min_card) || is.null(min_card))
      stop('min_card attribute must be specified')
    
    #compute maximum length
    maxLen <- max(sapply(Y0, nrow))
    full_data <- sapply(Y0, mtfd:::padding, maxLen,simplify = FALSE)
    
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
                          unsigned int portion_len){
  
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
#pragma omp parallel for
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
#pragma omp parallel for collapse(2)
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
    
    window_data_list <- createWindow(full_data,portion_len)
    window_data <- window_data_list[[1]]
    rownames(window_data) <- window_data_list[[2]]

    # set it as a matrix and omit rows with NAs
    window_data <- na.omit(window_data)
    
    # compute numerosity 
    numerosity <- rep(dim(Y0)[2] - portion_len + 1,times = dim(Y0)[1])
    removed_rows <- setdiff(window_data_list[[2]],rownames(window_data))
    curve_indices <- as.numeric(gsub("^(\\d+)_.*", "\\1", removed_rows))
    tab <- table(curve_indices)
    numerosity[as.numeric(names(tab))] <- numerosity[as.numeric(names(tab))] - tab
    
    rm(window_data_list)
    rm(removed_rows)
    rm(tab)
    
  cppFunction('
  Rcpp::NumericMatrix createDistance(Rcpp::NumericMatrix& windowData,
                                     const arma::vec& numerosity){
  
  int outrows = windowData.nrow();
  int outcols = windowData.ncol();
  bool isSingle = (numerosity.size() == 1);
  
  double outcols_inv = 1.0 / static_cast<double>(outcols);
  const arma::mat windowDataRef(windowData.begin(),outrows,outcols,false,true);
  const arma::colvec& vsum = arma::sum(windowDataRef,1) * outcols_inv;
  Rcpp::NumericMatrix result(outrows,outrows);
  arma::mat scoreData(result.begin(),outrows,outrows,false,true);
  
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
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
      scoreData(i,arma::span(start,end)) += M;
    }
    
    // Multiple curves
  } else {
    int until_here = 0;
    int numerosity_size = numerosity.size();
    for(int j = 0; j <numerosity_size;++j)
    {
      int num = numerosity[j] - 1;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
      for(int i = until_here; i < until_here + numerosity[j]-1;++i)
      {
        arma::uword start = i+1;
        arma::uword end = std::min<int>(std::min<int>(i+overlap ,i + (num--)),outrows);
        scoreData(arma::span(start,end),i) += M;
      }
      until_here += numerosity[j];
    }
  }
  return Rcpp::wrap(result);
}',depends = "RcppArmadillo")
# step 2: compute fMRS-based dissimilarity matrix
    D_fmsr <- createDistance(window_data,numerosity)
    rownames(D_fmsr) <- colnames(D_fmsr) <- rownames(window_data)
    ## STEP 3 -----
    # step 3: get the sub-trees (tree_s)
    minidend  <- mtfd:::get_minidend(as.dist(D_fmsr))
    # step 3: identify seeds and corresponding families

    all_paths   <- mtfd:::get_path_complete(minidend, window_data, min_card = min_card,worker_number = worker_number)
    all_paths   <- all_paths[!(lapply(all_paths, is.null) %>% unlist())] 
    
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

    ## STEP 4 -----
    # STEP 4: post-processing and rearranging results using different criteria
    
    ### CRITERION: adjusted fMSR ----
    best_order <- vec_of_scores %>% order()
    vec_of_scores_ordered <- vec_of_scores[best_order]
    # ordered list_of_recommendations
    list_of_recommendations_ordered <- list_of_recommendations[best_order]
 
    #for every recommended motif, compute all its accolites
    all_accolites <- lapply(list_of_recommendations_ordered,
                            function(x){
                              lapply(x, mtfd:::get_accolites, window_data, portion_len, FALSE) %>% unlist()
                            })
    
    #Starting from the top, we compare each motif to those with higher rank. If all portions of
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
    resFunBi = list(window_data = window_data,
                    list_of_recommendations_ordered = list_of_recommendations_ordered,
                    vec_of_scores_ordered = vec_of_scores_ordered)
    }
    else
    {
      resFunBi <- readRDS(paste0(name,'/resFunBi.rds'))
      window_data <- resFunBi$window_data
      list_of_recommendations_ordered <-  resFunBi$list_of_recommendations_ordered
      vec_of_scores_ordered <- resFunBi$vec_of_scores_ordered
    }

    # CRITERIUM: variance ----
    # We can order the results by variance too
    if(criterium == "Variance")
    {
      motif_var <- lapply(list_of_recommendations_ordered, 
                          function(x){
                            temp <- window_data[x,]
                            temp_mean <- temp %>% colMeans()
                            ((temp - temp_mean)^2 %>% sum())/(nrow(temp))
                          }) %>% unlist()
      
      var_order <- motif_var %>% order(decreasing = TRUE)
      list_of_recommendations_ordered <- list_of_recommendations_ordered[var_order]
      vec_of_scores_ordered <- vec_of_scores_ordered[var_order]
    }
    # save results 
    saveRDS(resFunBi,file = paste0(name,'/resFunBi.rds')) 
    # Creating some plot ----
    ## Plot the data and highlight the motif occurrences in red -----
    if(plot)
    {
      pdf(paste0(name,"/plot_",criterium,".pdf"),width=10,height=5)
      layout(matrix(1:2,ncol=2,byrow=TRUE),widths=c(8.5,1))
      for(q in 1:min(length(list_of_recommendations_ordered),cut_off,na.rm=TRUE)){
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
      for(q in 1:min(length(list_of_recommendations_ordered),cut_off,na.rm = TRUE)){
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
    return(list(list_of_recommendations_ordered = list_of_recommendations_ordered,
                vec_of_scores_ordered = vec_of_scores_ordered))
  }
}