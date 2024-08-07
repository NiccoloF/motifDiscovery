#' @title .find_min_diss
#'
#' @description Find shift warping minimizing dissimilarity between multidimensional curves (dimension=d).
#'
#' @param y list of two elements y0=y(x), y1=y'(x), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x), matrices with d columns.
#' @param alpha weight coefficient between d0.L2 and d1.L2.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param c_k minimum length of supp(y_shifted) and supp(v) intersection.
#' @return Shift warping and dissimilarity
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
.find_min_diss <- function(y,v,alpha,w,c_k,d,use0,use1,transform_y=FALSE,transform_v=FALSE){
  # Find shift warping minimizing dissimilarity between multidimensional
  # curves (dimension=d).
  # Return shift and dissimilarity.
  # y: list of two elements y0=y(x), y1=y'(x), matrices with d columns.
  # v: list of two elements v0=v(x), v1=v'(x), matrices with d columns.
  # alpha: weight coefficient between d0.L2 and d1.L2.
  # w: weights for the dissimilarity index in the
  # different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
  # transform_y: if TRUE, y is normalized to [0,1] before applying the distance.
  # transform_v: if TRUE, v is normalized to [0,1] before applying the distance.
  
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
