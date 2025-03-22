#' @title Find Minimum Dissimilarity
#' 
#' @description
#' Finds the shift warping that minimizes dissimilarity between multidimensional curves. 
#' The function operates on curves represented by two sets of data, each containing multiple dimensions.
#' 
#' @param y A list containing two matrices:
#'        - \code{y0}: The first set of curve values \(y(x)\).
#'        - \code{y1}: The second set of curve values \(y'(x)\).
#'        Each matrix should have \(d\) columns corresponding to the dimensions.
#' @param v A list containing two matrices:
#'        - \code{v0}: The first set of curve values \(v(x)\).
#'        - \code{v1}: The second set of curve values \(v'(x)\).
#'        Each matrix should have \(d\) columns corresponding to the dimensions.
#' @param alpha A numeric weight coefficient that balances the contributions of 
#'               the L2 norms of the two curve sets. 
#' @param w A numeric vector of weights for the dissimilarity index across different dimensions.
#'          All weights must be positive (\(w > 0\)).
#' @param c_k An integer specifying the minimum length of the intersection of the supports 
#'             of the shifted \(y\) and \(v\).
#' @param d An integer indicating the dimensionality of the curves.
#' @param use0 A logical value indicating whether to use the first component of the curves (i.e., \(y0\) and \(v0\)).
#' @param use1 A logical value indicating whether to use the second component of the curves (i.e., \(y1\) and \(v1\)).
#' @param transform_y A logical value indicating whether to normalize \(y\) to the range \([0,1]\) 
#'                    before calculating the distance.
#' @param transform_v A logical value indicating whether to normalize \(v\) to the range \([0,1]\) 
#'                    before calculating the distance.
#' 
#' @return A numeric vector containing:
#'         - The optimal shift that minimizes the dissimilarity.
#'         - The minimum dissimilarity value found.
#' 
#' @details
#' This function computes the shift warping between the provided multidimensional curves by 
#' examining various shifts and calculating the corresponding dissimilarity. The user can control 
#' which components of the curves to include in the calculation and whether to normalize the data.
#' 
#' The function returns both the optimal shift and the minimal dissimilarity, which can be used 
#' to assess the similarity between the two sets of curves under the specified constraints.
#' 
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
.find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1,transform_y=FALSE,transform_v=FALSE){
  v_dom=.domain(v,use0)
  v=.select_domain(v,v_dom,use0,use1)
  v_len=length(v_dom)
  y_len=unlist(lapply(y,nrow))[1]
  s_rep=(1-(v_len-c_k)):(y_len-v_len+1+(v_len-c_k))
  y_rep=lapply(s_rep,
               function(i){
                 index=i-1+seq_len(v_len)
                 y_rep_i=list(y0=NULL,y1=NULL)
                 if(use0){
                   y_rep_i$y0=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[1]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 if(use1){
                   y_rep_i$y1=rbind(matrix(NA,nrow=sum(index<=0),ncol=d),
                                    as.matrix(y[[2]][index[(index>0)&(index<=y_len)],]),
                                    matrix(NA,nrow=sum(index>y_len),ncol=d))
                 }
                 y_rep_i=.select_domain(y_rep_i,v_dom,use0,use1)
                 return(y_rep_i)
               })
  length_inter=unlist(lapply(y_rep,
                             function(y_rep_i){
                               if(use0)
                                 return(sum((!is.na(y_rep_i$y0[,1]))))
                               return(sum((!is.na(y_rep_i$y1[,1]))))
                             }))
  valid=length_inter>=c_k
  if(sum(valid)==0){
    valid[length_inter==max(length_inter)]=TRUE
  }
  s_rep=s_rep[valid]
  y_rep=y_rep[valid]
  d_rep=unlist(mapply(.diss_d0_d1_L2,y_rep,MoreArgs=list(v,w,alpha,transform_y,transform_v)))
  return(c(s_rep[which.min(d_rep)],min(d_rep)))
}
