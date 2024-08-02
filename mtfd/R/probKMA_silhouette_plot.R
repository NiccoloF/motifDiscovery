#' @title probKMA_silhouette
#'
#' @description Compute the adapted silhouette index on the results of probKMA.
#'
#' @param probKMA_results output of probKMA function (with return_options=TRUE).
#' @return A list containing:
#' @return \item{silhouette}{ vector of silhouette indices}
#' @return \item{motifs}{ vector of motifs numbers}
#' @return \item{curves}{ vector of curves numbers with motifs}
#' @return \item{silhouette_average}{ vector of average silhouette index for each cluster}
#' @author Marzia Angela Cremona  & Francesca Chiaromonte
probKMA_silhouette_plot <- function(silhouette_results,K,plot = TRUE){
  # Plot the adapted silhouette index on the results of probKMA.

  ### plot silhouette ########################################################################################
    if(plot) {
    silhouette = silhouette_results[[1]]
    Y_motifs = silhouette_results[[2]]
    curves_in_motifs = silhouette_results[[3]]
    silhouette_average = silhouette_results[[4]]
    curves_in_motifs_number = silhouette_results[[5]]
    
    
    n=length(silhouette)
    sil=rev(silhouette)
    y=barplot(sil,space=c(0,rev(diff(Y_motifs))),xlab='Silhouette index',names='',
              xlim=c(min(0,min(sil)-0.05),1.2),horiz=TRUE,las=1,mgp=c(2.5,1,0),col='lightblue')
    sapply(seq(0, 1, by = 0.2),function(h){abline(v = h, col = "firebrick3", lty = 2)})
    text(ifelse(sil>=0,sil+0.03,sil-0.03),ifelse(sil>0,y,y+0.2),labels=paste0("c",rev(unlist(curves_in_motifs))),cex=0.5)
    title(main='Silhouette plot',sub=paste("Average silhouette width:",round(mean(silhouette_average),2)),adj = 0.5)
    title(sub=paste("Average silhouette width:",round(mean(silhouette_average),2)),adj = 0.5,col.sub='red')
    mtext(paste("#curves =",n),adj=0,side = 3, line = 1,font=1,col=rgb(0,0.35,0))
    mtext(substitute("#motifs:"~K,list(K=K)),side = 3, line = 2, adj = 0,col=rgb(0,0.35,0))
    mtext(expression(paste(motif," | ",curve[j]," | avg ",s[i])),side = 3, line = 1.5, adj = 1,font=1,col=rgb(0,0.35,0))
    y=rev(y)
    for(k in seq_len(K)){
      y_k=mean(y[Y_motifs==k])
      text(1.1,y_k,paste0(k, " | ", curves_in_motifs_number[k], " | ", format(silhouette_average[k], digits = 1, nsmall = 2)),xpd=TRUE,srt=90,col='firebrick3')
    }
  }
  
  return(list(silhouette=silhouette_results[[1]],motifs=silhouette_results[[2]],curves=unlist(silhouette_results[[3]]),
              silhouette_average=silhouette_results[[4]]))
}