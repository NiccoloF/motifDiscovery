#' @title .find_diss
#'
#' @description Find dissimilarity between multidimensional curves (dimension=d), without alignment unless
#' their lengths are different. To be used by probKMA_silhouette fucntion.
#'
#' @param y list of two elements y0=y(x), y1=y'(x), matrices with d columns.
#' @param v list of two elements v0=v(x), v1=v'(x), matrices with d columns.
#' @param alpha weight coefficient between d0.L2 and d1.L2.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param aligned if TRUE, curves are already aligned. Else, the shortest curve is aligned inside the longest.
#' @return shift and dissimilarity
#' @author Marzia Angela Cremona & Francesca Chiaromonte
#' @export
.find_diss <- function(y,v,alpha,w,aligned,d,use0,use1){
  v_dom=.domain(v,use0)
  v=.select_domain(v,v_dom,use0,use1)
  v_len=length(v_dom)
  y_len=unlist(lapply(y,nrow))[1]
  if(aligned){
    s_rep=1
  }else{
    s_rep=1:(y_len-v_len+1)
  }
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
  d_rep=unlist(mapply(.diss_d0_d1_L2,y_rep,MoreArgs=list(v,w,alpha)))
  return(c(s_rep[which.min(d_rep)],min(d_rep)))
}
