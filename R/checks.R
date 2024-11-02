#' @title Initial Checks for ProbKMA
#' 
#' @description
#' This function performs various input checks on the parameters provided by the user to ensure they are valid 
#' for running the ProbKMA algorithm. It verifies the structure and content of the input data, including 
#' curves, derivatives, initial membership probabilities, shift matrices, and various parameters.
#' 
#' @param Y0 A list of matrices or vectors representing the curves.
#' @param Y1 A list of matrices or vectors representing the derivatives of the curves.
#' @param P0 A numeric matrix representing initial membership probabilities. Rows correspond to curves, and columns correspond to clusters.
#' @param S0 A numeric matrix representing the initial shift matrix.
#' @param params A list containing various parameters for ProbKMA, including:
#' \itemize{
#'    \item \code{standardize} Logical indicating whether to standardize the curves.
#'    \item \code{K} Number of motifs.
#'    \item \code{c} Minimum motif length.
#'    \item \code{c_max} Maximum motif length.
#'    \item \code{iter_max} Maximum number of iterations.
#'    \item \code{quantile} Quantile value for stopping criterion.
#'    \item \code{alpha} A numeric value related to the dissimilarity measure.
#'    \item \code{w} Weights used in the algorithm.
#'    \item \code{stopCriterion} Stopping criterion for the algorithm.
#'    \item \code{m} Weighting exponent.
#'    \item \code{tol} Tolerance level for stopping criteria.
#'    \item \code{iter4elong} Maximum iterations for elongation.
#'    \item \code{tol4elong} Tolerance for elongation.
#'    \item \code{max_elong} Maximum elongation allowed.
#'    \item \code{trials_elong} Number of trials for elongation.
#'    \item \code{deltaJK_elong} Threshold for elongation.
#'    \item \code{max_gap} Maximum gap allowed.
#'    \item \code{iter4clean} Number of iterations for cleaning.
#'    \item \code{tol4clean} Tolerance for cleaning.
#'    \item \code{quantile4clean} Quantile for cleaning.
#'    \item \code{return_options} Options for returning results.
#'    \item \code{seed} Seed for random number generation.
#'    \item \code{exe_print} Boolean to control printing of execution messages.
#'    \item \code{set_seed} Boolean to control whether to set a random seed.
#'    \item{\code{transformed}} A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.
#'    \item \code{n_threads} Number of threads for parallel processing.
#'  }
#' @param diss A character string indicating the type of dissimilarity measure to be used. 
#' Possible values are: \code{'d0_L2'}, \code{'d1_L2'}, \code{'d0_d1_L2'}.
#' @param V_init A list containing initial values for the clusters. If provided, it must match the expected structure based on \code{K}.
#' 
#' @return A list containing:
#' \item{FuncData}{A list of processed curves and derivatives after performing the checks.}
#' \item{Parameters}{A list of validated parameters ready for use in initializing the ProbKMA object.}
#' 
#' @export
initialChecks <- function(Y0,Y1,P0,S0,params,diss,V_init){
  ### Unpacking #############################################################################
  standardize = params$standardize
  K = params$K
  c = params$c
  c_max = params$c_max
  iter_max = params$iter_max
  quantile = params$quantile
  alpha = params$alpha
  w = params$w
  stop_criterion = params$stopCriterion
  m = params$m
  tol = params$tol
  iter4elong = params$iter4elong
  tol4elong = params$tol4elong
  max_elong = params$max_elong
  trials_elong = params$trials_elong
  deltaJk_elong = params$deltaJK_elong
  max_gap = params$max_gap
  iter4clean = params$iter4clean
  tol4clean = params$tol4clean
  quantile4clean = params$quantile4clean
  return_options = params$return_options
  seed = params$seed
  exe_print = params$exe_print
  set_seed = params$set_seed
  transformed = params$transformed
  n_threads = params$n_threads
  
  if(set_seed){
    set.seed(seed)
  }
  
  Y0_f <- function(y_i)
  {
    return(y_i[[1]])
  }
  Y1_f <- function(y_i)
  {
    return(y_i[[2]])
  }
  
  ### check input ####################################################################################
  start=proc.time()
  # check dissimilarity
  if(!(diss %in% c('d0_L2','d1_L2','d0_d1_L2')))
    stop('Dissimilarity not valid. Possible choices are: \'d0_L2\', \'d1_L2\' and \'d0_d1_L2\'.')
  # check Y0
  if(missing(Y0))
    stop('Y0 must be specified.')
  if(!is.null(Y0))
    if(!is.list(Y0))
      stop('Y0 should be a list of vectors or matrices.')
  if((FALSE %in% lapply(Y0,is.matrix))&&(FALSE %in% lapply(Y0,is.vector)))
    stop('Y0 should be a list of vectors or matrices.')
  N=length(Y0) # number of curves
  #if(N<5)
  #  stop('More curves y_i(x) needed.')
  Y0=lapply(Y0,as.matrix)
  d=unique(unlist(lapply(Y0,ncol),use.names=FALSE)) # dimension of curves
  if(length(d)>1)
    stop('Curves have not all the same dimension.')
  Y0_NA=lapply(Y0,function(y,d) sum(rowSums(!is.na(y)) %in% c(0,d))==nrow(y),d)
  if(FALSE %in% Y0_NA){
    warning('y_j(x)=NA only in some dimensions, for some x. Putting NA in all dimensions.')
    Y0=lapply(Y0,
              function(y,d){
                y[rowSums(!is.na(y))!=d,]=NA
                return(y)},
              d)
  }
  rm(Y0_NA)
  if(standardize){
    Y0_tot=Reduce(rbind,Y0)
    Y0_mean=colMeans(Y0_tot,na.rm=TRUE)
    Y0_sd=apply(Y0_tot,2,sd,na.rm=TRUE)
    rm(Y0_tot)
    Y0=lapply(Y0,function(y0,m,s) t((t(y0)-m)/s*100),m=Y0_mean,s=Y0_sd)
  }
  if(diss=='d0_d1_L2'){
    # check required input
    if(!is.numeric(alpha)){
      warning('alpha not valid. Setting alpha=0.5.')
      alpha=0.5
    }
    if((alpha<0)|(alpha>1)){
      warning('alpha not valid. Setting alpha=0.5.')
      alpha=0.5
    } else if(alpha==0){
      diss='d0_L2'
    } else if(alpha==1){
      diss='d1_L2'
    }
  }
  if(diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(Y0,function(y0) list(y0=y0,y1=NULL))
  }else{
    # check Y1
    if(is.null(Y1))
      stop(paste0('Y1 must be specified with diss=\'',diss,'\'.'))
    if(!is.list(Y1))
      stop('Y1 should be a list of vectors or matrices.')
    if((FALSE %in% lapply(Y1,is.matrix))&&(FALSE %in% lapply(Y1,is.vector)))
      stop('Y1 should be a list of vectors or matrices.')
    if(N!=length(Y1))
      stop('The number of curves in Y0 and derivatives Y1 should be the same.')
    Y1=lapply(Y1,as.matrix)
    d1=unique(unlist(lapply(Y1,ncol),use.names=FALSE)) # dimension of derivatives
    if(length(d1)>1)
      stop('Derivatives have not all the same dimension.')
    if(d!=d1)
      stop('The dimension of curves in Y0 and derivatives Y1 should be the same')
    Y1_NA=lapply(Y1,function(y,d) sum(rowSums(!is.na(y)) %in% c(0,d))==nrow(y),d)
    if(FALSE %in% Y1_NA){
      warning('y\'_j(x)=NA only in some dimensions, for some x. Putting NA in all dimensions.')
      Y1=lapply(Y1,
                function(y,d){
                  y[rowSums(!is.na(y))!=d,]=NA
                  return(y)},
                d)
    }
    rm(Y1_NA)
    Y0_NA=lapply(Y0,function(y) rowSums(is.na(y))!=0)
    Y1_NA=lapply(Y1,function(y) rowSums(is.na(y))!=0)
    diff_NA=mapply(function(y0,y1) y0!=y1,Y0_NA,Y1_NA)
    index_diff_NA=which(unlist(lapply(diff_NA,sum))!=0)
    if(length(index_diff_NA)>0){
      warning('y_j(x) and y\'_j(x) are not both defined, for some x. Putting NA in that case.')
      same_NA=mapply(function(y0,y1,diff_NA){
        y0[diff_NA]=NA
        y1[diff_NA]=NA
        return(list(y0,y1))
      },Y0[index_diff_NA],Y1[index_diff_NA],diff_NA[index_diff_NA])
      Y0[index_diff_NA]=same_NA[1,]
      Y1[index_diff_NA]=same_NA[2,]
    }
    rm(Y0_NA,Y1_NA,diff_NA)
    if(standardize){
      Y1=lapply(Y1,function(y1,s) t(t(y1)/s*100),s=Y0_sd)
    }
    if(diss=='d1_L2'){
      alpha=1
      use0=FALSE
      use1=TRUE
      Y=lapply(Y1,function(y1) list(y0=NULL,y1=y1))
    }
    if(diss=='d0_d1_L2'){
      use0=TRUE
      use1=TRUE
      Y=mapply(function(y0,y1) list(y0=y0,y1=y1),Y0,Y1,SIMPLIFY=FALSE)
    }
  }
  # check required input
  if(missing(K))
    stop('K must be specified')
  if(missing(c))
    stop('c must be specified')
  # check K
  if(length(K)!=1)
    stop('Number of motifs K not valid.')
  if(K%%1!=0)
    stop('Number of motifs K should be an integer.')
  if(K<1)
    stop('Number of motifs K should be at least 1.')
  # check c
  if(!(length(c) %in% c(1,K)))
    stop('Minimum motif length c not valid.')
  if(sum(c%%1!=0))
    stop('Minimum motif lengths should be integers.')
  if(sum(c<=1))
    stop('Minimum motif lengths should be at least 2.')
  
  #check V_init
  if(!is.null(V_init)){
    if(length(V_init)!=K){
      V_init=NULL
      warning('The length of the list does not represent the number of clusters. The random
         initialization will be used.')
    }
    for(i in 1:length(V_init)){
      if((nrow(V_init[[i]][[1]])!=c)&(nrow(V_init[[i]][[2]])!=c)){ 
        #use nrow instead of length so it will work in a d-dimensional case as well.
        V_init=NULL
        warning('The length of the list does not represent the number of clusters. The random
         initialization will be used.')  
      }
    }
    
  }
  
  if(!is.null(V_init)){
    V_init <- list(lapply(V_init,Y0_f), lapply(V_init,Y1_f))
  }
  
  c=rep(c,length.out=K)
  # find all intervals contained in supp(y_i) and their lengths
  Y_intervals=lapply(Y0,
                     function(y,d){
                       y_not_NA=(rowSums(!is.na(y))==d) # find supp(y)
                       intervals=which((y_not_NA[2:length(y_not_NA)]-y_not_NA[1:(length(y_not_NA)-1)])==1)+1
                       if(y_not_NA[1])
                         intervals=c(1,intervals)
                       intervals_lengths=rle(y_not_NA)$lengths[rle(y_not_NA)$values==TRUE]
                       intervals_all=data.frame(start=intervals,
                                                end=intervals+intervals_lengths-1,
                                                length=intervals_lengths)
                       return(intervals_all)
                     },d)
  # check minimum motif length compatibility
  min_length=unlist(lapply(Y_intervals,
                           function(y_intervals,c){
                             return(sum(y_intervals$length>=c))},
                           max(c)))
  if(0 %in% min_length)
    stop('Minimum motifs length is incompatible with supplied curves. Choose a smaller minimum motifs length.')
  # check c_max
  if(!(length(c_max) %in% c(1,K)))
    stop('Maximum motif length c_max not valid.')
  c_max=rep(c_max,length.out=K)
  for(k in 1:K){
    if(c_max[k]!=Inf){
      if(c_max[k]%%1!=0)
        stop('Maximum motif lengths should be integers.')
      if(c_max[k]<=1)
        stop('Maximum motif length should be at least 2.')
      # check maximum motif length compatibility
      max_length=unlist(lapply(Y_intervals,
                               function(y_intervals,c_max){
                                 return(sum(floor((1+max_gap)*y_intervals$length)>=c_max))},
                               c_max[k]))
      if(0 %in% max_length){
        warning('Maximum motif length is incompatible with supplied curves. Choosing default maximum motif length for motif ',k,'.')
        c_max[k]=Inf
      }
    }
    if(c_max[k]==Inf){
      c_max[k]=floor((1+max_gap)*min(unlist(lapply(Y_intervals,function(y_intervals) max(y_intervals$length)))))
    }
  }
  # check P0
  if(!is.null(P0)){
    if(sum(dim(P0)!=c(N,K))!=0){
      warning('Membership matrix dimensions not valid. Choosing random initial membership matrix.')
      P0=NULL
    }else{
      if(sum((P0<0)|(P0>1))){
        warning('Memberships should be non-negative and <=1. Choosing random initial membership matrix.')
        P0=NULL
      }else{
        if(sum(rowSums(P0)!=1)){
          warning('Memberships of each curve should sum to 1. Choosing random initial membership matrix.')
          P0=NULL
        }else{
          if(sum(colSums(P0)==0)){
            warning('Sum of memberships of each cluster should be positive. Choosing random initial membership matrix.')
            P0=NULL
          }
        }
      }
    }
  }
  # check S0
  if(!is.null(S0)){
    if(sum(dim(S0)!=c(N,K))!=0){
      warning('Shift warping matrix dimensions not valid. Choosing random initial shift warping matrix.')
      S0=NULL
    }else{
      Y_segments=mapply(function(y,S_i,c){
        return(mapply(function(s,c) NA %in% y[s+seq_len(c)-1],S_i,c))},
        Y0,lapply(seq_len(N),function(i) S0[i,]),MoreArgs=list(c),SIMPLIFY=TRUE)
      if(sum(Y_segments)){
        warning('Shift warping matrix not valid. Choosing random initial shift warping matrix.')
        S0=NULL
      }
    }
  }
  # check weigths
  if(!is.vector(w)|!is.numeric(w))
    stop('Weights w not valid.')
  if(!(length(w) %in% c(1,d)))
    stop('Weights w not valid.')
  if(TRUE %in% (w<=0))
    stop('Weights w not valid.')
  w=rep(w,length.out=d)
  # check weighting exponent
  if(!is.numeric(m)|(length(m)>1))
    stop('Weighting exponent m not valid.')
  if(m<=1)
    stop('Weighting exponent m should be >1')
  # check maximum number of iterations
  if(!is.numeric(iter_max)|(length(iter_max)!=1))
    stop('Maximum number of iterations iter_max not valid.')
  if(((iter_max%%1)!=0)|(iter_max<=0))
    stop('Maximum number of iterations iter_max not valid.')
  # check stop criterion
  if(!(stop_criterion %in% c('max','mean','quantile')))
    stop('Stop criterion not valid. Possible choices are: \'max\', \'mean\', \'quantile\'.')
  if(stop_criterion=='quantile')
    if((!is.numeric(quantile))||(length(quantile)!=1)||(quantile<0)||(quantile>1))
      stop('quantile should be a number in [0,1].')
  # check tolerance
  if((!is.numeric(tol))||(length(tol)!=1)||(tol<=0))
    stop('Tolerance should be a positive number.')
  # check elongation input
  if((!is.numeric(iter4elong))||(length(iter4elong)!=1))
    stop('iter4elong not valid.')
  if(((iter4elong%%1)!=0)|(iter4elong<=0))
    stop('iter4elong not valid.')
  if((!is.numeric(tol4elong))||(length(tol4elong)!=1)||(tol4elong<=0))
    stop('tol4elong should be a positive number.')
  if((!is.numeric(max_elong))||(length(max_elong)!=1)||(max_elong<=0))
    stop('max_elong should be a positive number.')
  if((!is.numeric(trials_elong))||(length(trials_elong)!=1))
    stop('trials_elong not valid.')
  if(((trials_elong%%1)!=0)|(trials_elong<=0))
    stop('trials_elong not valid.')
  if((!is.numeric(deltaJk_elong))||(length(deltaJk_elong)!=1)||(deltaJk_elong<=0))
    stop('deltaJk_elong should be a positive number.')
  if((!is.numeric(max_gap))||(length(max_gap)!=1)||(max_gap<0))
    stop('max_gap should be a non-negative number.')
  # check cleaning input
  if((!is.numeric(iter4clean))||(length(iter4clean)!=1))
    stop('iter4clean not valid.')
  if(((iter4clean%%1)!=0)|(iter4clean<=0))
    stop('iter4clean not valid.')
  if((!is.numeric(tol4clean))||(length(tol4clean)!=1)||(tol4clean<=0))
    stop('tol4clean should be a positive number.')
  if((!is.numeric(quantile4clean))||(length(quantile4clean)!=1)||(N*K*quantile4clean<K)||(quantile4clean>1))
    stop('quantile4clean not valid.')
 
  ### initialize #############################################################################################
  start=proc.time()
  if(is.null(P0)){
    # create random membership matrix, with N rows and k columns
    P0=matrix(c(runif(N*(K-1)),rep(1,N)),nrow=N,ncol=K)
    P0=as.matrix(apply(P0,1,sort))
    if(K>1)
      P0=cbind(P0[1,],Reduce('cbind',lapply(2:K,function(k) P0[k,]-P0[k-1,])))
  }
  colnames(P0)=NULL
  if(is.null(S0)){
    # create shift warping matrix, with N rows and k columns
    S0=matrix(unlist(lapply(Y_intervals,
                            function(y_intervals,c,K){
                              s0=rep(NA,K)
                              for(k in 1:K){
                                y_intervals_k=y_intervals[y_intervals$length>=c[k],]
                                y_starts=unlist(apply(y_intervals_k,1,
                                                      function(y_interval)
                                                        return(y_interval[1]:(y_interval[2]-c[k]+1))),
                                                use.names=FALSE)
                                s0[k]=sample(y_starts,1)
                              }
                              return(s0)
                            },c,K)),
              nrow=N,ncol=K,byrow=TRUE)
  }
  
  Y <- list("Y0" = lapply(Y,Y0_f),"Y1" = lapply(Y,Y1_f))

  return(list("FuncData" = list("Y"=Y,"P0"=P0,"S0"=S0,"V_init"=V_init),
              "Parameters" = list("standardize"=standardize,"K"=K,"c"=c,"c_max"=c_max,
                                  "iter_max"=iter_max,"quantile"=quantile,
                                  "stopCriterion"=stop_criterion,"m"=m,"w"=w,
                                  "alpha"=alpha,"tol"=tol,
                                  "iter4elong"=iter4elong,"tol4elong"=tol4elong,
                                  "max_elong"=max_elong,"trials_elong"=trials_elong,
                                  "deltaJK_elong"=deltaJk_elong,"max_gap"=max_gap,
                                  "iter4clean"=iter4clean,"tol4clean"=tol4clean,
                                  "quantile4clean"=quantile4clean,
                                  "seed"=seed,
                                  "return_options"=return_options,
                                  "exe_print"= exe_print,
                                  "set_seed" = set_seed,
                                  "transformed" = transformed,
                                  "n_threads" = n_threads) ) )
}


