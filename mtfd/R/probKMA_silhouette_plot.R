#' @title probKMA_silhouette
#'
#' @description Compute the adapted silhouette index on the results of probKMA.
#'
#' @param probKMA_results output of probKMA function (with return_options=TRUE).
#' @param align if TRUE, try all possible alignements between pieces of curves (corresponding to the
#' same or to different motifs).
#' @param plot if TRUE, the silhouette plot is drawn.
#' @return A list containing:
#' @return \item{silhouette}{ vector of silhouette indices}
#' @return \item{motifs}{ vector of motifs numbers}
#' @return \item{curves}{ vector of curves numbers with motifs}
#' @return \item{silhouette_average}{ vector of average silhouette index for each cluster}
#' @author Marzia Angela Cremona  & Francesca Chiaromonte
#' @export
probKMA_silhouette_plot <- function(results, align=FALSE, plot=TRUE){
  # Compute the adapted silhouette index on the results of probKMA.
  # probKMA_results: output of probKMA function (with return_options=TRUE).
  # align: if TRUE, try all possible alignements between pieces of curves (corresponding to the same or to different motifs).
  # plot: if TRUE, the silhouette plot is drawn.
  
  # check silhouette_results exists
  if(length(results) != 2){
    stop('probKMA executed without probKMA_silhouette function')
  }
  
  probKMA_results = results[[1]]
  silhouette_results = results[[2]]
  
  ### plot silhouette ########################################################################################
  if(plot){
    K=ncol(probKMA_results$P_clean)
    n=length(silhouette_results[[1]])
    sil=rev(silhouette_results[[1]])
    y=barplot(sil,space=c(0,rev(diff(silhouette_results[[2]]))),xlab='Silhouette index',names='',
              xlim=c(min(0,min(sil)-0.05),1.2),horiz=TRUE,las=1,mgp=c(2.5,1,0),col='gray')
    text(ifelse(sil>=0,sil+0.03,sil-0.03),ifelse(sil>0,y,y+0.2),labels=rev(unlist(silhouette_results[[3]])),cex=0.5)
    title(main='Silhouette plot')
    title(sub=paste("Average silhouette width:",round(mean(silhouette_results[[4]]),2)),adj=0)
    mtext(paste("n =",n),adj=0)
    mtext(substitute(K~~"motifs",list(K=K)),adj=1)
    mtext(expression(paste(j," : ",n[j]," | avg ",s[i])),adj=1.04,line=-1.2)
    y=rev(y)
    for(k in seq_len(K)){
      y_k=mean(y[silhouette_results[[2]]==k])
      text(1,y_k,
           paste0(k,": ",silhouette_results[[5]][k]," | ",format(silhouette_results[[4]][k],digits=1,nsmall=2)),xpd=NA,adj=0.1)
    }
  }
  
  return(list(silhouette=silhouette_results[[1]],motifs=silhouette_results[[2]],curves=unlist(silhouette_results[[3]]),
              silhouette_average=silhouette_results[[4]]))
}