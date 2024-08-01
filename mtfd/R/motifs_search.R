#' @title motifs_search
#'
#' @description Find occurrences of the candidate motifs in the curves and sort 
#' them according to their frequencies and radius,
#' In each group (as defined by cutting the dendrogram at high 2*Rall), 
#' we choose the motif
#' with highest frequency and lower mean dissimilarity (the one ranking best in both dimensions),
#' Additional motifs can be chosen in a group, if their lengths differ enough from the length of the first motif chosen,
#' A candidate motif matches a piece of curve if their dissimilarity is less than the corresponding R_m.
#'
#' @param cluster_candidate_motifs_results output of cluster_candidate_motifs function.
#' @param R_all global radius, used to cut the dendrogram (requiring groups to be more than 2*Rall apart).
#' @param R_m vector with group-specific radii, used to find motif occurrences. The length of the vector
#' must match the number of clusters obtained cutting the dendrogram at height 2*Rall. If NULL, Rm is
#' determined in each group (based on distances between motifs of the same group and all curves).
#' @param use_real_occurrences if TRUE, find occurrences for all candidate motifs and uses real frequency
#' and mean dissimilarity to choose motifs in groups (more accurate, but time consuming). Otherwise,
#' uses approximate frequency and mean dissimilarity (default).
#' @param length_diff minimum difference in length among motifs of the same group, required in ordered
#' to keep more than one motif, in percentage of the most frequent motif.
#' @param worker_number number of CPU cores to be used for parallelization (default number of CPU cores).
#' If worker_number=1, the function is run sequentially.
#' @return A list containing Y0 and Y1 output of cluster_candidate_motifs function and:
#' @return \item{V0}{ list of founds motifs}
#' @return \item{V1}{ list of derived of founds motifs}
#' @return \item{V_length}{ vector of real lengths of founds motifs}
#' @return \item{V_occurrences}{ vector of real occurrences of founds motifs}
#' @return \item{V_frequencies}{ vector of real frequencies of founds motifs}
#' @return \item{V_mean_diss}{ vector of real average dissimilarities between the founds motifs}
#' @return \item{R_motifs}{ vector of radii of founds motifs}
#' @author Marzia Angela Cremona & Francesca Chiaromonte
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
  library(parallel)
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
        pred_knn=knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=cluster_candidate_motifs_results$k_knn,prob=TRUE)
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
        pred_knn=knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=cluster_candidate_motifs_results$k_knn,prob=TRUE)
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