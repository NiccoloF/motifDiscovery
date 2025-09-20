#' @title Plot Silhouette Index from probKMA Results
#'
#' @description This function generates a bar plot displaying the adapted silhouette index based on the results from the `probKMA` algorithm. It visually represents the quality of the motifs identified by illustrating the average silhouette width for each motif. Additionally, it provides relevant information about the number of curves associated with each motif.
#'
#' @param silhouette_results A list containing the results from the silhouette analysis, including:
#'   - Silhouette indices for each motif.
#'   - Motif identifiers.
#'   - Curves associated with each motif.
#'   - Average silhouette widths.
#'   - Number of curves in each motif.
#' @param K An integer representing the number of motifs identified by the `probKMA` algorithm.
#' @param plot A logical value indicating whether to generate the plot. Default is `TRUE`.
#'
#' @return A list containing the following elements:
#' \item{silhouette}{A vector of silhouette indices for the motifs.}
#' \item{motifs}{A vector of motif identifiers.}
#' \item{curves}{A vector containing all curves associated with the motifs.}
#' \item{silhouette_average}{A vector of average silhouette widths for each motif.}
#'
#' @import graphics
#' @import grDevices
#' @export
probKMA_silhouette_plot <- function(silhouette_results,K,plot = TRUE){
  # Plot the adapted silhouette index on the results of probKMA.
  oldpar <- par(no.readonly = TRUE)
  # Ensure the original settings are restored when the function exits
  on.exit(par(oldpar))
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