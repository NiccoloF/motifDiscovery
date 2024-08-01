#' @title .diss_d0_d1_L2
#'
#' @description Dissimilarity index for multidimensional curves (dimension=d).
#' Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
#'
#' @param y list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param alpha weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
#' @export
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
.diss_d0_d1_L2 <- function(y,v,w,alpha,transform_y=FALSE,transform_v=FALSE){
  # Dissimilarity index for multidimensional curves (dimension=d).
  # Sobolev type distance with normalization on common support: (1-alpha)*d0.L2+alpha*d1.L2.
  # y: list of two elements y0=y(x), y1=y'(x) for x in dom(v), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x) for x in dom(v), matrices with d columns.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # alpha: weight coefficient between d0.L2 and d1.L2 (alpha=0 means d0.L2, alpha=1 means d1.L2).
  # transform_y: if TRUE, y is normalized to [0,1] before applying the distance.
  # transform_v: if TRUE, v is normalized to [0,1] before applying the distance.
  
  y_norm=y
  v_norm=v
  
  if(transform_y){
    y0_min = apply(y[[1]], 2, min, na.rm = TRUE)
    y0_max = apply(y[[1]], 2, max, na.rm = TRUE)
    y0_diff = y0_max - y0_min
    y0_const = unlist(lapply(y0_diff, function(diff) all.equal(diff, 0) == TRUE))
    y_norm[[1]] = t( (t(y[[1]]) - y0_min) / y0_diff )
    y_norm[[2]] = t( t(y[[2]]) / y0_diff )
    y_norm[[1]][,y0_const]=0.5
    y_norm[[2]][,y0_const] = 0
  } 
  if(transform_v){
    v0_min=apply(v[[1]],2,min,na.rm=TRUE)
    v0_max=apply(v[[1]],2,max,na.rm=TRUE)
    v0_diff=v0_max-v0_min
    v0_const=unlist(lapply(v0_diff,function(diff) all.equal(diff,0)==TRUE))
    v_norm[[1]]=t( (t(v[[1]]) - v0_min)/v0_diff)
    v_norm[[2]]=t( t(v[[2]])/v0_diff)
    v_norm[[1]][,v0_const]=0.5
    v_norm[[2]][,v0_const]=0
  }
  .diss_L2 <- function(y,v,w){
    # Dissimilarity index for multidimensional curves (dimension=d).
    # L2 distance with normalization on common support.
    # y=y(x), v=v(x) for x in dom(v), matrices with d columns.
    # w: weights for the dissimilarity index in the different dimensions (w>0).
    
    sum(colSums((y-v)^2,na.rm=TRUE)/(colSums(!is.na(y)))*w)/ncol(y)
    # NB: divide for the length of the interval, not for the squared length!
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