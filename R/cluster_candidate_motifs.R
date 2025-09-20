#' @title Cluster Candidate Motifs
#'
#' @description This function clusters candidate motifs based on their distances and computes group-specific radii for motif clusters. It utilizes K-nearest neighbors (KNN) for determining a global radius and evaluates overlaps among motifs. The function supports parallel computation for efficiency.
#'
#' @param filter_candidate_motifs_results A list containing results from filtering candidate motifs, including various components like `Y0`, `Y1`, `V0_clean`, `V1_clean`, `D_clean`, and more, which are essential for the clustering process.
#' @param motif_overlap A numeric value representing the minimum proportion of overlap required between motifs to be considered similar (default is 0.6).
#' @param k_knn An integer specifying the number of nearest neighbors to consider when determining the global radius (default is 3).
#' @param votes_knn_Rall A numeric value indicating the threshold for KNN voting when determining the global radius (default is 0.5).
#' @param votes_knn_Rm A numeric value indicating the threshold for KNN voting when determining group-specific radii (default is 0.5).
#' @param worker_number An optional integer specifying the number of parallel workers to use. If NULL, it defaults to the number of available cores minus one.
#'
#' @return A list containing:
#' - `VV_D`: Matrix of distances between motifs.
#' - `VV_S`: Matrix of shifts between motifs.
#' - `k_knn`: The value of K used in KNN.
#' - `votes_knn_Rall`: Voting threshold for the global radius.
#' - `R_all`: The global radius determined from the clustering process.
#' - `hclust_res`: Result of hierarchical clustering (if applicable).
#' - `votes_knn_Rm`: Voting threshold for group-specific radius.
#' - `R_m`: Vector of group-specific radii for each cluster.
#' - All components from the input `filter_candidate_motifs_results`.
#'
#' @details This function performs the following steps:
#' 1. Sets up parallel jobs based on the specified `worker_number`.
#' 2. Prepares input data based on the type of distance measure used.
#' 3. Computes distances between motifs.
#' 4. Determines a global radius (`R_all`) using KNN classification.
#' 5. Clusters motifs and determines group-specific radii (`R_m`) for each cluster.
#'
#' @importFrom combinat combn
#' @importFrom stats as.dist quantile density na.omit as.dendrogram
#' @importFrom utils modifyList
#' @import ggplot2
#' @export
cluster_candidate_motifs <- function(filter_candidate_motifs_results,motif_overlap=0.6,
                                     k_knn=3,votes_knn_Rall=0.5,votes_knn_Rm=0.5,worker_number=NULL){
 
  ### set parallel jobs ###################################################################################
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
  if(worker_number>1){
    cl_search=parallel::makeCluster(worker_number,timeout=60*60*24*30)
    on.exit(parallel::stopCluster(cl_search))
  }else{
    cl_search=NULL
  }
  
  ### prepare input data ##################################################################################
  if(filter_candidate_motifs_results$diss=='d0_L2'){
    alpha=0
    use0=TRUE
    use1=FALSE
    Y=lapply(filter_candidate_motifs_results$Y0,function(y0) list(y0=y0,y1=NULL))
    V=lapply(filter_candidate_motifs_results$V0_clean,function(v0) list(v0=v0,v1=NULL))
  }
  if(filter_candidate_motifs_results$diss=='d1_L2'){
    alpha=1
    use0=FALSE
    use1=TRUE
    Y=lapply(filter_candidate_motifs_results$Y1,function(y1) list(y0=NULL,y1=y1))
    V=lapply(filter_candidate_motifs_results$V1_clean,function(v1) list(v0=NULL,v1=v1))
  }
  if(filter_candidate_motifs_results$diss=='d0_d1_L2'){
    alpha=filter_candidate_motifs_results$alpha
    use0=TRUE
    use1=TRUE
    Y=mapply(function(y0,y1) list(y0=y0,y1=y1),filter_candidate_motifs_results$Y0,filter_candidate_motifs_results$Y1,SIMPLIFY=FALSE)
    V=mapply(function(v0,v1) list(v0=v0,v1=v1),filter_candidate_motifs_results$V0_clean,filter_candidate_motifs_results$V1_clean,SIMPLIFY=FALSE)
  }
  w=filter_candidate_motifs_results$w
  max_gap=filter_candidate_motifs_results$max_gap
  d=ncol(filter_candidate_motifs_results$Y0[[1]])
  N=nrow(filter_candidate_motifs_results$D)
  K=ncol(filter_candidate_motifs_results$D)
  V_dom=lapply(filter_candidate_motifs_results$V0,function(v) rowSums(!is.na(v))!=0)
  V_length=unlist(lapply(V_dom,length))
  
  
  ### compute distances between motifs ####################################################################
  if(length(V)==1){
    VV_D=as.matrix(0)
    VV_S=as.matrix(1)
  }else{
    VV=combn(V,2,simplify=FALSE)
    VV=array(unlist(VV,recursive=FALSE),dim=c(2,length(VV)))
    VV_lengths=as.matrix(combn(V_length,2))
    VV_motif_overlap=floor(apply(VV_lengths,2,min)*motif_overlap)
    SD=mapply(.find_min_diss,VV[1,],VV[2,],VV_motif_overlap,
              MoreArgs=list(alpha=alpha,w=w,d=d,use0=use0,use1=use1),SIMPLIFY=TRUE)
    VV_D=matrix(0,nrow=length(V),ncol=length(V))
    VV_D[lower.tri(VV_D)]=SD[2,]
    VV_D=VV_D+t(VV_D) # matrix of distances
    VV_S=matrix(0,nrow=length(V),ncol=length(V))
    VV_S[lower.tri(VV_S)]=SD[1,]
    VV_S=VV_S+t(VV_S)+diag(1,nrow=length(V)) # matrix of shifts
  }
  
  ### determine a global radius Rall ######################################################################
  dataR=data.frame(P=as.vector(filter_candidate_motifs_results$P_clean),D=as.vector(filter_candidate_motifs_results$D_clean))
  if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
    R_all=max(dataR$D[dataR$P==1])
  }else{
    D_new=seq(0,max(filter_candidate_motifs_results$D_clean),length.out=10000)
    pred_knn=class::knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=k_knn,prob=TRUE)
    R_all=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=votes_knn_Rall)[1]]
  }
  
  ### cluster motifs and determine a group-specific radius Rm ##############################################
  if(length(V)==1){
    hclust_res=NULL
    R_m=R_all
  }else{
    VV_dist=as.dist(VV_D)
    hclust_res=hclust(VV_dist,method='average') # hierarchical clustering based on motif-motif distances
    V_hclust=cutree(hclust_res,h=2*R_all) # cut at high 2*R_all
    n_hclust=max(V_hclust)
    
    R_m=rep(NA,n_hclust)
    for(i_hclust in seq_len(n_hclust)){
      index_i=which(V_hclust==i_hclust) # motifs in cluster i
      V_length_i=V_length[index_i]
      V_P_i=as.matrix(filter_candidate_motifs_results$P_clean[,index_i])
      V_D_i=as.matrix(filter_candidate_motifs_results$D_clean[,index_i])
      
      dataR=data.frame(P=as.vector(V_P_i),D=as.vector(V_D_i))
      if(max(dataR$D[dataR$P==1])<=min(dataR$D[dataR$P==0])){
        R_m[i_hclust]=max(dataR$D[dataR$P==1])
      }else{
        D_new=seq(0,max(V_D_i),length.out=10000)
        pred_knn=class::knn(train=dataR$D,test=as.matrix(D_new),cl=dataR$P,k=k_knn,prob=TRUE)
        R_m[i_hclust]=D_new[which(ifelse(pred_knn==1,attributes(pred_knn)$prob,1-attributes(pred_knn)$prob)<=votes_knn_Rm)[1]]
      }
    }
  }
  
  ### output ###############################################################################################
  return(c(list(VV_D=VV_D,VV_S=VV_S,k_knn=k_knn,votes_knn_Rall=votes_knn_Rall,R_all=R_all,hclust_res=hclust_res,votes_knn_Rm=votes_knn_Rm,R_m=R_m),filter_candidate_motifs_results)) 
}