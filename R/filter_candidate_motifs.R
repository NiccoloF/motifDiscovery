#' @title Filter Candidate Motifs
#' 
#' @description
#' Filters the candidate motifs based on a specified threshold for the average silhouette index and 
#' a threshold for the size of the curves in the motif. This function is useful for refining the set 
#' of candidate motifs by removing those that do not meet the defined criteria.
#' 
#' @param find_candidate_motifs_results Output from the \code{find_candidate_motif} function, which contains the results of the motif discovery process.
#' @param sil_threshold A numeric value representing the threshold for the average silhouette index. 
#'                      Values should be between -1 and 1.
#' @param size_threshold An integer representing the threshold for the size of the motif, defined as the number of curves in the cluster. 
#'                       Should be at least 1.
#' @param K A vector containing the numbers of motifs that must be considered. 
#'          This should be a subset of \code{find_candidate_motifs_results$K}.
#' @param c A vector with the minimum motif lengths that must be considered. 
#'          This should be a subset of \code{find_candidate_motifs_results$c}.
#' 
#' @return A list containing the filtered candidate motifs and some ProbKMA options:
#' \item{V0_clean}{A vector of candidate motifs after filtering.}
#' \item{V1_clean}{A vector of derived candidate motifs after filtering.}
#' \item{D_clean}{A matrix of distances of candidate motifs from the curves after filtering.}
#' \item{P_clean}{A matrix of probabilities of membership of the candidate motifs after filtering.}
#' \item{Y0}{The original input data for the first curve.}
#' \item{Y1}{The original input data for the derivative of the curve.}
#' \item{diss}{The dissimilarity matrix used in motif discovery.}
#' \item{alpha}{The weight parameter for the Sobolev distance.}
#' \item{w}{Weights for the dissimilarity index across different dimensions.}
#' \item{max_gap}{Maximum allowable gap between curves.}
#' 
#' @details
#' This function checks the provided thresholds for silhouette indices and motif sizes, 
#' ensuring they are valid before filtering the motifs. If any parameters are invalid, 
#' default values from the \code{find_candidate_motifs_results} are used.
#' 
#' The function then loads the necessary data files for each combination of motif numbers and sizes, 
#' filtering them using the silhouette criteria. The resulting motifs are sorted by their length 
#' in descending order, and the filtered results are returned in a structured list.
#' 
#' @export
filter_candidate_motifs <- function(find_candidate_motifs_results,sil_threshold=0.5,size_threshold=2,
                                    K=find_candidate_motifs_results$K,c=find_candidate_motifs_results$c){
  
  ### check input #############################################################################################
  # check sil_threshold
  if(length(sil_threshold)!=1)
    stop('sil_threshold not valid.')
  if((sil_threshold<(-1))|(sil_threshold>1))
    stop('sil_threshold should be a number bewteen -1 and 1.')
  # check size_threshold
  if(length(size_threshold)!=1)
    stop('size_threshold not valid.')
  if(size_threshold%%1!=0)
    stop('size_threshold should be an integer.')
  if(size_threshold<1)
    stop('size_threshold should be at least 1.')
  # check K
  if(TRUE %in% !(K %in% find_candidate_motifs_results$K)){
    warning('K is not a subset of find_candidate_motifs_results$K. Using default K.')
    K=find_candidate_motifs_results$K
  }
  # check c
  if(TRUE %in% !(c %in% find_candidate_motifs_results$c)){
    warning('c is not a subset of find_candidate_motifs_results$c. Using default c.')
    c=find_candidate_motifs_results$c
  }
  
  ### filter motifs ###########################################################################################  
  n_init=find_candidate_motifs_results$n_init
  name=find_candidate_motifs_results$name
  motifs=lapply(K,
                function(K){
                  lapply(c,
                         function(c){
                           lapply(seq_len(n_init),
                                  function(i){
                                    load(paste0(name,"K",K,"_c",c,'/random',i,'.RData'))
                                    motifs=probKMA_silhouette_filter(probKMA_results,silhouette,sil_threshold,size_threshold)
                                    return(motifs)
                                  })
                         })
                })
  motifs=unlist(unlist(motifs,recursive=FALSE),recursive=FALSE)
  if(is.null(unlist(motifs)))
    stop("No motif present after filtering. Please re-run the function with less stringent parameters.")
  V0_clean=unlist(lapply(motifs,function(motifs) motifs$V0_clean),recursive=FALSE)
  V_clean_length=unlist(lapply(V0_clean,length))
  index=order(V_clean_length,decreasing=TRUE) # order from the longest to the shortest
  V0_clean=V0_clean[index]
  V1_clean=unlist(lapply(motifs,function(motifs) motifs$V1_clean),recursive=FALSE)[index]
  D_clean=Reduce(cbind,lapply(motifs,function(motifs) motifs$D_clean))[,index]
  P_clean=Reduce(cbind,lapply(motifs,function(motifs) motifs$P_clean))[,index]
  c=Reduce(c,lapply(motifs,function(motifs) motifs$c))[index]
  K=Reduce(c,lapply(motifs,function(motifs) motifs$K))[index]
  
  ### output ##################################################################################################
  load(paste0(name,"K",K[1],"_c",c[1],'/random',1,'.RData'))
  return(list(V0_clean=V0_clean,V1_clean=V1_clean,D_clean=as.matrix(D_clean),P_clean=as.matrix(P_clean),c=c,K=K,
              Y0=probKMA_results$Y0,Y1=probKMA_results$Y1,
              diss=probKMA_results$diss,alpha=probKMA_results$alpha,w=probKMA_results$w,max_gap=probKMA_results$max_gap))
}