#' @title Find Occurrences of a Motif
#' 
#' @description
#' Finds occurrences of a specified motif in a set of curves, where the dissimilarity is lower than a specified threshold \( R \). 
#' The function compares the motif against curves represented in multiple dimensions and returns the details of matching occurrences.
#' 
#' @param v A list containing two matrices:
#'          - \code{v0}: The first set of motif values.
#'          - \code{v1}: The second set of motif values.
#'          Each matrix should have \( d \) columns corresponding to the dimensions.
#' @param Y A list of \( N \) lists, each containing two matrices:
#'          - \code{Y0}: The first set of curve values.
#'          - \code{Y1}: The second set of curve values.
#'          Each matrix should have \( d \) columns corresponding to the dimensions.
#' @param R A numeric value representing the maximum allowed dissimilarity.
#' @param alpha A numeric value that serves as a weight coefficient between the L2 norms 
#'               of the two sets of motifs when using the dissimilarity function \code{diss_d0_d1_L2}.
#' @param w A numeric vector of weights for the dissimilarity index across different dimensions.
#'          All weights must be positive (\( w > 0 \)).
#' @param c_k An integer specifying the minimum length of the intersection of the supports 
#'             of the shifted motif and the curves.
#' @param use0 A logical value indicating whether to use the first component of the curves (i.e., \( Y0 \) and \( v0 \)).
#' @param use1 A logical value indicating whether to use the second component of the curves (i.e., \( Y1 \) and \( v1 \)).
#' @param transformed A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.
#' 
#' @return A matrix with three columns:
#'         - \code{curve}: The ID of the curve where the motif was found.
#'         - \code{shift}: The optimal shift at which the motif occurs.
#'         - \code{diss}: The dissimilarity value associated with the match.
#'         If no occurrences are found, an empty matrix is returned.
#' 
#' @details
#' The function systematically checks each curve in the provided list against the specified motif 
#' to identify positions where the dissimilarity does not exceed the defined threshold \( R \). 
#' It handles multidimensional curves and can selectively consider different components of the motifs and curves.
#' 
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
.find_occurrences <- function(v,Y,R,alpha,w,c_k,use0,use1,transformed=FALSE){
  v_dom=.domain(v,use0)
  v_len=length(v_dom)
  SD_motif=lapply(Y,
                  function(y){
                    y_len=unlist(lapply(y,nrow))[1]
                    s_rep=seq_len(y_len-v_len+1)
                    y_rep=lapply(s_rep,
                                 function(i){
                                   y_rep=list(y0=NULL,y1=NULL)
                                   if(use0)
                                     y_rep$y0=as.matrix(y$y0[i-1+seq_len(v_len),])
                                   if(use1)
                                     y_rep$y1=as.matrix(y$y1[i-1+seq_len(v_len),])
                                   y_rep=.select_domain(y_rep,v_dom,use0,use1)
                                 })
                    valid=unlist(lapply(lapply(y_rep,.domain,use0),sum))>=c_k
                    s_rep=s_rep[valid]
                    y_rep=y_rep[valid]
                    d_rep=unlist(lapply(y_rep,.diss_d0_d1_L2,.select_domain(v,v_dom,use0,use1),w,alpha,transformed))
                    d_rep_R=c(FALSE,d_rep<=R,FALSE)
                    diff_d_rep_R=diff(d_rep_R)
                    start=which(diff_d_rep_R==1)
                    end=which(diff_d_rep_R==(-1))-1
                    SD_motif=mapply(function(start,end){
                      index=(start:end)[which.min(d_rep[start:end])]
                      return(c(s_rep[index],d_rep[index]))
                    },start,end)
                    return(SD_motif)
                  })
  if(!is.null(unlist(SD_motif))){
    v_occurrences=cbind(rep(seq_along(SD_motif),unlist(lapply(SD_motif,length))/2),
                        matrix(unlist(SD_motif),ncol=2,byrow=TRUE))
    row.names(v_occurrences)=NULL
    colnames(v_occurrences)=c('curve','shift','diss')
  }else{
    v_occurrences=c()
  }
  return(v_occurrences)
}
