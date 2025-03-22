#' @title Motif Search in Curves
#'
#' @description
#' The `motifs_search` function identifies and ranks motifs within a set of curves based on their 
#' frequencies and dissimilarity measures. It processes candidate motifs clustered from hierarchical 
#' clustering results, selects optimal motifs within each cluster, and determines their occurrences 
#' in the original curves. The function supports parallel processing to enhance computational efficiency 
#' and offers flexibility in handling different dissimilarity metrics and motif selection criteria.
#'
#' @param cluster_candidate_motifs_results A list containing the output from the `cluster_candidate_motifs` 
#'   function. This list must include elements such as:
#'   \describe{
#'     \item{Y0}{A list of matrices representing the original curves.}
#'     \item{Y1}{A list of matrices representing the derivatives of the curves (if applicable).}
#'     \item{V0_clean}{A list of candidate motifs derived from `Y0`.}
#'     \item{V1_clean}{A list of candidate motifs derived from `Y1` (if applicable).}
#'     \item{D_clean}{A matrix of dissimilarity measures between motifs and curves.}
#'     \item{P_clean}{A matrix indicating positive matches (e.g., presence of motifs in curves).}
#'     \item{hclust_res}{A hierarchical clustering object obtained from `hclust`.}
#'     \item{R_all}{A numeric value representing the global radius used for dendrogram cutting.}
#'     \item{w}{A numeric vector of weights for the dissimilarity index across different dimensions.}
#'     \item{transformed}{A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.}
#'     \item{max_gap}{A numeric value defining the maximum allowed gap in distances for cluster separation.}
#'     \item{k_knn}{An integer specifying the number of neighbors for K-Nearest Neighbors classification.}
#'     \item{votes_knn_Rm}{A numeric value defining the probability threshold for KNN-based radius determination.}
#'     \item{c}{A numeric vector specifying the minimum number of overlapping elements required for motif validation.}
#'   }
#' 
#' @param R_all A numeric value representing the global radius used to cut the dendrogram, ensuring that 
#'   clusters are at least twice this radius apart. This parameter defines the grouping of motifs into clusters.
#'
#' @param R_m A numeric vector containing group-specific radii used to identify motif occurrences within each cluster. 
#'   The length of this vector must match the number of clusters obtained by cutting the dendrogram at a height of `2 * R_all`. 
#'   If `NULL`, the function automatically determines `R_m` for each group based on the distances between motifs 
#'   within the same cluster and all curves.
#'
#' @param different_R_m_finding A logical value indicating whether to use a different radius (`R_m_finding`) 
#'   for finding motif occurrences compared to the initial radius (`R_m`). If `TRUE`, `R_m_finding` is used; 
#'   otherwise, `R_m` is employed. This allows for separate tuning of motif occurrence detection.
#'
#' @param R_m_finding A numeric vector containing group-specific radii used specifically for finding motif occurrences 
#'   when `different_R_m_finding` is set to `TRUE`. The length of this vector must match the number of clusters obtained 
#'   by cutting the dendrogram at a height of `2 * R_all`. If `NULL`, `R_m_finding` is determined automatically for 
#'   each group based on distances between motifs within the same cluster and all curves.
#'
#' @param use_real_occurrences A logical value indicating whether to compute real occurrences of candidate motifs 
#'   within the curves. If `TRUE`, the function calculates actual frequencies and mean dissimilarities for motif 
#'   selection, providing more accurate results at the cost of increased computation time. If `FALSE`, it uses 
#'   approximate frequencies and mean dissimilarities for faster execution. Defaults to `FALSE`.
#'
#' @param length_diff A numeric value specifying the minimum percentage difference in length required among 
#'   motifs within the same group to retain multiple motifs. This parameter ensures diversity in motif selection 
#'   by preventing motifs of similar lengths from being selected simultaneously. It is defined as a percentage 
#'   relative to the length of the most frequent motif. Defaults to `Inf`, meaning no additional motifs are selected 
#'   based on length differences.
#'
#' @param worker_number An integer indicating the number of CPU cores to utilize for parallel processing. 
#'   By default, the function uses one less than the total number of available cores (`detectCores() - 1`). 
#'   Setting `worker_number = 1` forces the function to run sequentially without parallelization. If `NULL`, 
#'   the function automatically determines the optimal number of workers based on the system's available cores.
#'
#' @return A list containing:
#'   \describe{
#'     \item{V0}{A list of selected motifs derived from `Y0`.}
#'     \item{V1}{A list of selected motifs derived from `Y1` (if applicable).}
#'     \item{V_length}{A numeric vector representing the real lengths of the selected motifs.}
#'     \item{V_occurrences}{A list detailing the occurrences of each selected motif within the curves.}
#'     \item{V_frequencies}{A numeric vector indicating the real frequencies of each selected motif.}
#'     \item{V_mean_diss}{A numeric vector representing the average dissimilarity of each selected motif.}
#'     \item{Y0}{A list of matrices corresponding to the original curves, as provided in `cluster_candidate_motifs_results`.}
#'     \item{Y1}{A list of matrices corresponding to the derivatives of the curves (if applicable), as provided in `cluster_candidate_motifs_results`.}
#'     \item{R_motifs}{A numeric vector containing the radii associated with each selected motif.}
#'   }
#'
#' @details
#' The `motifs_search` function operates through the following steps:
#' \enumerate{
#'   \item **Parallelization Setup**: Determines the number of worker cores to use based on `worker_number`. If `worker_number > 1`, 
#'     it initializes a cluster for parallel processing.
#'   \item **Input Preparation**: Depending on the dissimilarity metric (`d0_L2`, `d1_L2`, or `d0_d1_L2`), it prepares the data structures 
#'     `Y` and `V` for processing.
#'   \item **Dendrogram Cutting**: Cuts the hierarchical clustering dendrogram at a height of `2 * R_all` to define clusters of motifs.
#'   \item **Radius Determination**: If `R_m` or `R_m_finding` is not provided, the function calculates these radii for each cluster 
#'     based on motif distances and K-Nearest Neighbors (KNN) classification.
#'   \item **Candidate Motif Selection**: Depending on `use_real_occurrences`, the function either computes real occurrences and uses 
#'     actual frequencies and mean dissimilarities to select motifs, or it uses approximate measures for faster processing.
#'   \item **Motif Filtering**: Within each cluster, motifs are ranked based on their frequency and mean dissimilarity. Additional motifs 
#'     can be selected if their lengths differ sufficiently from the most frequent motif, as defined by `length_diff`.
#'   \item **Output Compilation**: The selected motifs and their associated properties are compiled into a comprehensive list for further analysis or visualization.
#' }
#'
#' @import parallel
#' @importFrom dplyr %>%
#' @importFrom data.table as.data.table setnames
#' @importFrom class knn
#' @importFrom fastcluster hclust
#' @importFrom Rcpp evalCpp
#' @export
motifs_search <- function(cluster_candidate_motifs_results,
                          R_all=cluster_candidate_motifs_results$R_all,R_m=NULL,different_R_m_finding=FALSE,R_m_finding=NULL,
                          use_real_occurrences=FALSE,length_diff=Inf,worker_number=NULL){
  # Find occurrences of the candidate motifs in the curves and sort them according to their frequencies and radius.
  # In each group (as defined by cutting the dendrogram at high 2*Rall), we choose the motif
  # with highest frequency and lower mean dissimilarity (the one ranking best in both dimensions).
  # Additional motifs can be chosen in a group, if their lengths differ enough from the length of the first motif chosen.
  # A candidate motif matches a piece of curve if their dissimilarity is less than the corresponding R_m.
  # cluster_candidate_motifs_results: output of cluster_candidate_motifs function.
  # R_all: global radius, used to cut the dendrogram (requiring groups to be more than 2*Rall apart).
  # R_m: vector with group-specific radii, used to find motif occurrences. The length of the vector must match the number of clusters
  #      obtained cutting the dendrogram at height 2*Rall. If NULL, Rm is determined in each group (based on distances between motifs of the same group and all curves).
  # different_R_m_finding: if TRUE, find the final occurrences using the R_m_finding radius instead of the R_m
  # R_m_finding: vector with group-specific radii, used to find motif occurrences if different_R_m_finding= TRUE. The length of the vector must match the number of clusters
  #      obtained cutting the dendrogram at height 2*Rall. If NULL, Rm for finding is determined in each group (based on distances between motifs of the same group and all curves).
  # use_real_occurrences: if TRUE, find occurrences for all candidate motifs and uses real frequency and mean dissimilarity to choose motifs
  #                       in groups (more accurate, but time consuming). Otherwise, uses approximate frequency and mean dissimilarity (default).
  # length_diff: minimum difference in length among motifs of the same group, required in ordered to keep more than one motif, in % of the most frequent motif.
  # worker_number: number of CPU cores to be used for parallelization (default number of CPU cores). If worker_number=1, the function is run sequentially.
  
  ### set parallel jobs ###################################################################################
  core_number <- detectCores()
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
  if(worker_number>1){
    cl_search=makeCluster(worker_number,timeout=60*60*24*30)
    on.exit(stopCluster(cl_search))
  }else{
    cl_search=NULL
  }
  
  ### prepare input data ##################################################################################
  if(cluster_candidate_motifs_results$diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(cluster_candidate_motifs_results$Y0,function(y0) list(y0=y0,y1=NULL))
    V=lapply(cluster_candidate_motifs_results$V0_clean,function(v0) list(v0=v0,v1=NULL))
  }
  if(cluster_candidate_motifs_results$diss=='d1_L2'){
    alpha=1
    use0=FALSE
    use1=TRUE
    Y=lapply(cluster_candidate_motifs_results$Y1,function(y1) list(y0=NULL,y1=y1))
    V=lapply(cluster_candidate_motifs_results$V1_clean,function(v1) list(v0=NULL,v1=v1))
  }
  if(cluster_candidate_motifs_results$diss=='d0_d1_L2'){
    alpha=cluster_candidate_motifs_results$alpha
    use0=TRUE
    use1=TRUE
    Y=mapply(function(y0,y1) list(y0=y0,y1=y1),cluster_candidate_motifs_results$Y0,cluster_candidate_motifs_results$Y1,SIMPLIFY=FALSE)
    V=mapply(function(v0,v1) list(v0=v0,v1=v1),cluster_candidate_motifs_results$V0_clean,cluster_candidate_motifs_results$V1_clean,SIMPLIFY=FALSE)
  }
  w=cluster_candidate_motifs_results$w
  transformed=cluster_candidate_motifs_results$transformed
  max_gap=cluster_candidate_motifs_results$max_gap
  d=ncol(cluster_candidate_motifs_results$Y0[[1]])
  N=nrow(cluster_candidate_motifs_results$D)
  K=ncol(cluster_candidate_motifs_results$D)
  V_dom=lapply(cluster_candidate_motifs_results$V0,function(v) rowSums(!is.na(v))!=0)
  V_length=unlist(lapply(V_dom,length))
  
  ### cut dendrogram and check group-specific radius Rm ##################################################
  V_hclust=cutree(cluster_candidate_motifs_results$hclust_res,h=2*R_all) # cut at high 2*R_all
  n_hclust=max(V_hclust)
  
  # re-compute or check group-specific radius Rm and Rm for finding
  if(is.null(R_m)){
    R_m=rep(NA,n_hclust)
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_length_i=V_length[index_i]
      V_P_i=as.matrix(cluster_candidate_motifs_results$P_clean[,index_i])
      V_D_i=as.matrix(cluster_candidate_motifs_results$D_clean[,index_i])
      
      dataR=data.frame(P=as.vector(V_P_i),D=as.vector(V_D_i))
      if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
        R_m[i_hclust]=max(dataR$D[dataR$P==1])
      }else{
        D_new=seq(0,max(V_D_i),length.out=10000)
        pred_knn=class::knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=cluster_candidate_motifs_results$k_knn,prob=TRUE)
        R_m[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=cluster_candidate_motifs_results$votes_knn_Rm)[1]]
      }
    }
  }
  if(length(R_m)!=n_hclust)
    stop(paste0('The length of the vector R_m must match the number of clusters: ',n_hclust))
  if(is.null(R_m_finding)){
    R_m_finding=rep(NA,n_hclust)
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_length_i=V_length[index_i]
      V_P_i=as.matrix(cluster_candidate_motifs_results$P_clean[,index_i])
      V_D_i=as.matrix(cluster_candidate_motifs_results$D_clean[,index_i])
      
      dataR=data.frame(P=as.vector(V_P_i),D=as.vector(V_D_i))
      if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
        R_m_finding[i_hclust]=max(dataR$D[dataR$P==1])
      }else{
        D_new=seq(0,max(V_D_i),length.out=10000)
        pred_knn=class::knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=cluster_candidate_motifs_results$k_knn,prob=TRUE)
        R_m_finding[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=cluster_candidate_motifs_results$votes_knn_Rm_finding)[1]]
      }
    }
  }
  if(length(R_m_finding)!=n_hclust)
    stop(paste0('The length of the vector R_m_finding must match the number of clusters: ',n_hclust))
  
  ### select candidate motifs and find occurrences #######################################################
  if(use_real_occurrences){
    ### find candidate motifs in the curves ####################################################################
    c_k=floor(V_length*(1-max_gap))
    c_k[c_k<cluster_candidate_motifs_results$c]=cluster_candidate_motifs_results$c[c_k<cluster_candidate_motifs_results$c]
    V_R_m=R_m[V_hclust]
    V_R_m_finding=R_m_finding[V_hclust]
    # find occurrences
    V_occurrences=.mapply_custom(cl_search,.find_occurrences,
                                 V,V_R_m,c_k,MoreArgs=list(Y=Y,alpha=alpha,w=w,use0,use1,transformed),SIMPLIFY=FALSE)
    not_null=which(!unlist(lapply(V_occurrences,is.null)))
    V=V[not_null]
    V_length=V_length[not_null]
    V_R_m=V_R_m[not_null]
    V_R_m_finding=V_R_m_finding[not_null]
    V_hclust=V_hclust[not_null]
    
    ### select candidate motifs in each group ##############################################################
    select=c()
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      
      V_frequencies_i=unlist(lapply(V_occurrences[index_i],nrow)) # real frequency
      V_mean_diss_i=unlist(lapply(V_occurrences[index_i],function(x) mean(x[,3]))) # real average distance (actual radius)
      # order based of frequency and average distance
      V_order_i=order(rank(-V_frequencies_i)+rank(V_mean_diss_i)) # sum of ranks in the two dimensions
      V_frequencies_i=V_frequencies_i[V_order_i]
      V_mean_diss_i=V_mean_diss_i[V_order_i]
      index_i_ordered=index_i[V_order_i]
      V_length_i=V_length[index_i_ordered]
      
      # select motifs to keep
      keep=rep_len(TRUE,length(index_i))
      for(i in seq_len(length(index_i))){
        if(keep[i]){
          select=c(select,index_i_ordered[i])
          keep[i]=FALSE
          # motifs with length different enough from length of selected motif
          length_diff_enough=keep&((V_length_i<(V_length_i[i]*(1-length_diff)))|(V_length_i>(V_length_i[i]*(1+length_diff))))
          keep[!length_diff_enough]=FALSE
        }
      }
    }
    V=V[select]
    V_length=V_length[select]
    V_R_m=V_R_m[select]
    V_R_m_finding=V_R_m_finding[select]
    
    if(!different_R_m_finding){
      V_occurrences=V_occurrences[select]
      V_frequencies=unlist(lapply(V_occurrences,nrow)) # real frequency
      V_mean_diss=unlist(lapply(V_occurrences,function(x) mean(x[,3]))) # real average distance (actual radius)
      # sort final motifs based of frequency and average distance
      V_order=order(rank(-V_frequencies)+rank(V_mean_diss)) # sum of ranks in the two dimensions
      V=V[V_order]
      V_occurrences=V_occurrences[V_order]
      V_length=V_length[V_order]
      V_R_m=V_R_m[V_order]
      V_R_m_finding=V_R_m_finding[V_order]
      V_frequencies=V_frequencies[V_order]
      V_mean_diss=V_mean_diss[V_order]
      
      index_final=not_null[select][V_order]
    }else{
      ### find candidate motifs in the curves ####################################################################
      c_k=floor(V_length*(1-max_gap))
      c_k[c_k<cluster_candidate_motifs_results$c[select]]=cluster_candidate_motifs_results$c[select][c_k<cluster_candidate_motifs_results$c[select]]
      # find occurrences
      V_occurrences=.mapply_custom(cl_search,.find_occurrences,
                                   V,V_R_m_finding,c_k,MoreArgs=list(Y=Y,alpha=alpha,w=w,use0,use1,transformed),SIMPLIFY=FALSE)
      not_null=which(!unlist(lapply(V_occurrences,is.null)))
      V=V[not_null]
      V_occurrences=V_occurrences[not_null]
      V_length=V_length[not_null]
      V_R_m=V_R_m[not_null]
      V_R_m_finding=V_R_m_finding[not_null]
      V_frequencies=unlist(lapply(V_occurrences,nrow)) # real frequency
      V_mean_diss=unlist(lapply(V_occurrences,function(x) mean(x[,3]))) # real average distance (actual radius)
      # sort final motifs based of frequency and average distance
      V_order=order(rank(-V_frequencies)+rank(V_mean_diss)) # sum of ranks in the two dimensions
      V=V[V_order]
      V_occurrences=V_occurrences[V_order]
      V_length=V_length[V_order]
      V_R_m=V_R_m[V_order]
      V_R_m_finding=V_R_m_finding[V_order]
      V_frequencies=V_frequencies[V_order]
      V_mean_diss=V_mean_diss[V_order]
      
      index_final=select[not_null][V_order]
    }
  }else{
    ### select candidate motifs in each group ##############################################################
    select=c()
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_D_i=as.matrix(cluster_candidate_motifs_results$D_clean[,index_i])
      
      V_frequencies_approx_i=colSums(V_D_i<=R_m[i_hclust]) # approximate frequency
      V_mean_diss_approx_i=apply(V_D_i,2,function(x) mean(x[x<=R_m[i_hclust]])) # approximate average distance (actual radius)
      # order based of frequency and average distance
      V_order_i=order(rank(-V_frequencies_approx_i)+rank(V_mean_diss_approx_i)) # sum of ranks in the two dimensions
      V_frequencies_approx_i=V_frequencies_approx_i[V_order_i]
      V_mean_diss_approx_i=V_mean_diss_approx_i[V_order_i]
      index_i_ordered=index_i[V_order_i]
      V_length_i=V_length[index_i_ordered]
      
      # select motifs to keep
      keep=rep_len(TRUE,length(index_i))
      for(i in seq_len(length(index_i))){
        if(keep[i]){
          select=c(select,index_i_ordered[i])
          keep[i]=FALSE
          # motifs with length different enough from length of selected motif
          length_diff_enough=keep&((V_length_i<(V_length_i[i]*(1-length_diff)))|(V_length_i>(V_length_i[i]*(1+length_diff))))
          keep[!length_diff_enough]=FALSE
        }
      }
    }
    V=V[select]
    V_length=V_length[select]
    V_R_m=R_m[V_hclust[select]]
    V_R_m_finding=R_m_finding[V_hclust[select]]
    
    ### find candidate motifs in the curves ####################################################################
    c_k=floor(V_length*(1-max_gap))
    c_k[c_k<cluster_candidate_motifs_results$c[select]]=cluster_candidate_motifs_results$c[select][c_k<cluster_candidate_motifs_results$c[select]]
    # find occurrences
    if(different_R_m_finding){
      V_occurrences=.mapply_custom(cl_search,.find_occurrences,
                                   V,V_R_m_finding,c_k,MoreArgs=list(Y=Y,alpha=alpha,w=w,use0,use1,transformed),SIMPLIFY=FALSE)
    }else{
      V_occurrences=.mapply_custom(cl_search,.find_occurrences,
                                   V,V_R_m,c_k,MoreArgs=list(Y=Y,alpha=alpha,w=w,use0,use1,transformed),SIMPLIFY=FALSE)
    }
    not_null=which(!unlist(lapply(V_occurrences,is.null)))
    V=V[not_null]
    V_occurrences=V_occurrences[not_null]
    V_length=V_length[not_null]
    V_R_m=V_R_m[not_null]
    V_R_m_finding=V_R_m_finding[not_null]
    V_frequencies=unlist(lapply(V_occurrences,nrow)) # real frequency
    V_mean_diss=unlist(lapply(V_occurrences,function(x) mean(x[,3]))) # real average distance (actual radius)
    # sort final motifs based of frequency and average distance
    V_order=order(rank(-V_frequencies)+rank(V_mean_diss)) # sum of ranks in the two dimensions
    V=V[V_order]
    V_occurrences=V_occurrences[V_order]
    V_length=V_length[V_order]
    V_R_m=V_R_m[V_order]
    V_R_m_finding=V_R_m_finding[V_order]
    V_frequencies=V_frequencies[V_order]
    V_mean_diss=V_mean_diss[V_order]
    
    index_final=select[not_null][V_order]
  }
  
  ### output ##################################################################################################
  if(cluster_candidate_motifs_results$diss=='d0_L2'){
    V0=lapply(V,function(v) v$v0)
    V1=cluster_candidate_motifs_results$V1_clean[index_final]
  }else if(cluster_candidate_motifs_results$diss=='d1_L2'){
    V0=cluster_candidate_motifs_results$V0_clean[index_final]
    V1=lapply(V,function(v) v$v1)
  }else{
    V0=lapply(V,function(v) v$v0)
    V1=lapply(V,function(v) v$v1)
  }
  if(different_R_m_finding){
    R_motifs=V_R_m_finding
  }else{
    R_motifs=V_R_m
  }
  return(list(V0=V0,V1=V1,
              V_length=V_length,V_occurrences=V_occurrences,V_frequencies=V_frequencies,V_mean_diss=V_mean_diss,
              Y0=cluster_candidate_motifs_results$Y0,Y1=cluster_candidate_motifs_results$Y1,R_motifs=R_motifs))
}