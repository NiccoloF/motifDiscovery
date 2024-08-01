#' @title find_occurrences
#'
#' @description Find occurrences of a motif in a set of curves (dimesion=d), with dissimilarity lower than R.
#'
#' @param v list of 2 elements, v0, v1, matrices with d columns.
#' @param Y list of N lists of two elements, Y0, Y1, matrices with d columns.
#' @param R maximum dissimilarity allowed.
#' @param alpha if diss_fun=diss_d0_d1_L2, weight coefficient between d0_L2 and d1_L2.
#' @param w weights for the dissimilarity index in the different dimensions (w>0).
#' @param c_k minimum length of supp(y_shifted) and supp(v) intersection.
#' @return curve id, shift and dissimilarity
#' @author Marzia Angela Cremona & Francesca Chiaromonte
.find_occurrences <- function(v,Y,R,alpha,w,c_k,use0,use1,transformed=FALSE){
  # Find occurrences of a motif in a set of curves (dimesion=d), with dissimilarity lower than R.
  # Return curve id, shift and dissimilarity.
  # v: list of 2 elements, v0, v1, matrices with d columns.
  # Y: list of N lists of two elements, Y0, Y1, matrices with d columns.
  # R: maximum dissimilarity allowed.
  # alpha: if diss_fun=diss_d0_d1_L2, weight coefficient between d0_L2 and d1_L2.
  # w: weights for the dissimilarity index in the different dimensions (w>0).
  # c_k: minimum length of supp(y_shifted) and supp(v) intersection.
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
