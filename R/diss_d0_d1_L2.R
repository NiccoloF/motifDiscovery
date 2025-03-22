#' @title Dissimilarity Index for Multidimensional Curves
#' 
#' @description
#' Computes a Sobolev-type dissimilarity index for multidimensional curves based on a weighted combination of L2 norms of the function and its derivative. 
#' The distance is normalized on a common support, allowing for a comparison between curves considering both their levels and variations.
#' 
#' @param y A list of two matrices: 
#'   \itemize{
#'     \item \code{y[[1]]}: Values of the curve at points in the domain (y(x)).
#'     \item \code{y[[2]]}: Values of the derivative of the curve at the same points (y'(x)).
#'   }
#'   Each matrix should have \code{d} columns, where \code{d} is the dimensionality of the curves.
#' 
#' @param v A list of two matrices: 
#'   \itemize{
#'     \item \code{v[[1]]}: Values of the reference curve at points in the domain (v(x)).
#'     \item \code{v[[2]]}: Values of the derivative of the reference curve at the same points (v'(x)).
#'   }
#'   Each matrix should have \code{d} columns, where \code{d} is the dimensionality of the curves.
#' 
#' @param w A numeric vector of weights for the dissimilarity index in different dimensions. Each weight should be greater than 0.
#' 
#' @param alpha A numeric value (between 0 and 1) that specifies the weight coefficient between the L2 norm of the function (d0.L2) and the L2 norm of the derivative (d1.L2):
#'   \itemize{
#'     \item \code{alpha = 0}: Only the levels (d0.L2) are considered.
#'     \item \code{alpha = 1}: Only the derivative information (d1.L2) is considered.
#'     \item Values between 0 and 1 provide a mixture of both.
#'   }
#' 
#' @param transform_y A logical value indicating whether to normalize \code{y} to the range [0, 1] before applying the distance computation. Default is \code{FALSE}.
#' 
#' @param transform_v A logical value indicating whether to normalize \code{v} to the range [0, 1] before applying the distance computation. Default is \code{FALSE}.
#' 
#' @return A numeric value representing the dissimilarity index between the curves defined by \code{y} and \code{v}.
#' 
#' @details
#' The dissimilarity index is calculated based on the following Sobolev-type distance:
#' \deqn{D = (1 - \alpha) \cdot d0.L2 + \alpha \cdot d1.L2}
#' where:
#' \itemize{
#'   \item \code{d0.L2}: L2 distance considering only the levels of the curves.
#'   \item \code{d1.L2}: L2 distance considering only the derivatives of the curves.
#' }
#' 
#' The function normalizes the inputs based on the specified flags to ensure that all features are comparable.
#' 
#' @export
.diss_d0_d1_L2 <- function(y,v,w,alpha,transform_y=FALSE,transform_v=FALSE){

  y_norm=y
  v_norm=v
  
  .diss_L2 <- function(y,v,w){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # L2 distance with normalization on common support.
    # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    
    sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y)
    # NB: divide for the length of the interval, not for the squared length!
  }
  
  if(transform_y){
    y0_min = apply(y[[1]], 2, min, na.rm = TRUE)
    y0_max = apply(y[[1]], 2, max, na.rm = TRUE)
    y0_diff = y0_max - y0_min
    y0_const = unlist(lapply(y0_diff, function(diff) all.equal(diff, 0) == TRUE))
    y_norm[[1]] = t( (t(y[[1]]) - y0_min) / y0_diff )
    y_norm[[1]][,y0_const]=0.5
    if(alpha != 0) {
      y_norm[[2]] = t( t(y[[2]]) / y0_diff )
      y_norm[[2]][,y0_const] = 0
    }
  } 
  if(transform_v){
    v0_min=apply(v[[1]],2,min,na.rm=TRUE)
    v0_max=apply(v[[1]],2,max,na.rm=TRUE)
    v0_diff=v0_max-v0_min
    v0_const=unlist(lapply(v0_diff,function(diff) all.equal(diff,0)==TRUE))
    v_norm[[1]]=t( (t(v[[1]]) - v0_min)/v0_diff)
    v_norm[[1]][,v0_const]=0.5
    if(alpha != 0) {
      v_norm[[2]]=t( t(v[[2]])/v0_diff)
      v_norm[[2]][,v0_const]=0
    }
  }
  
  if(alpha==0){#L2-like distance which focuses excusively on the levels
    return(.diss_L2(y_norm[[1]],v_norm[[1]],w))
  }else if(alpha==1){#L2-like pseudo distance, which uses only weak
    #derivative information
    return(.diss_L2(y_norm[[2]],v_norm[[2]],w))
  }else{#Sobolev-like distance that allows to highlight
    #more complex features of curve shapes,taking into account
    #both levels and variations
    return((1-alpha)*.diss_L2(y_norm[[1]],v_norm[[1]],w)+alpha*.diss_L2(y_norm[[2]],v_norm[[2]],w))
  }
}