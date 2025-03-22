#' @title Plot Motif Search Results
#'
#' @description
#' The `motifs_search_plot` function visualizes the results obtained from the `motifs_search` function.
#' It generates plots for detected motifs across multiple dimensions, displaying both the motifs and their
#' corresponding derivative curves (if available). Users can filter motifs based on frequency thresholds
#' and choose to display either all motifs or the top `n` motifs. Additionally, the function provides an
#' option to plot all underlying curves with colored motifs highlighted.
#'
#' @param motifs_search_results A list containing the output from the `motifs_search` function. This
#'   includes elements such as `V0`, `V1`, `V_frequencies`, `Y0`, `Y1`, `V_length`, `V_occurrences`,
#'   `V_mean_diss`, and `R_motifs`, which store information about the detected motifs and their properties.
#'
#' @param ylab A character string specifying the label for the y-axis in the plots. This label will be
#'   appended with the dimension number to create individual titles for each plot. Default is an empty
#'   string (`''`).
#'
#' @param freq_threshold An integer indicating the minimum frequency a motif must have to be included
#'   in the plots. Only motifs with a frequency equal to or greater than `freq_threshold` will be
#'   visualized. Default value is `5`.
#'
#' @param top_n Determines how many motifs to plot based on their frequency. If set to `'all'`, all
#'   motifs meeting the `freq_threshold` will be plotted. If an integer is provided, only the top
#'   `top_n` motifs with the highest frequencies will be displayed. Default is `'all'`.
#'
#' @param plot_curves A logical value indicating whether to plot all underlying curves with colored
#'   motifs highlighted. If `TRUE`, the function overlays motifs on the curves for better visualization.
#'   Default is `TRUE`.
#'
#' @param transformed A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.
#'
#' @return
#' The function does not return a value but generates plots visualizing the motifs and their occurrences
#' across different dimensions. It creates separate plots for each dimension and includes legends for
#' easy identification of motifs.
#'
#' @details
#' The `motifs_search_plot` function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters, ensuring that the frequency threshold and `top_n` are appropriate.
#'   \item Selects motifs that meet the frequency threshold and, if specified, limits the number of motifs
#'     to the top `n`.
#'   \item For each dimension, it plots the motif centers (`V0`) and their derivatives (`V1`, if available).
#'   \item If `plot_curves` is `TRUE`, it overlays the motifs on the original curves, highlighting them with
#'     distinct colors.
#'   \item Adds legends to the plots for clear identification of each motif.
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom data.table as.data.table setnames
#' @export
motifs_search_plot <- function(motifs_search_results,ylab='',freq_threshold=5,top_n='all',plot_curves=TRUE,transformed = FALSE){
  # Plot the results of motifs_search.
  # motifs_search_results: output of motifs_search function.
  # freq_threshold: plot only motifs with frequency at least equal to freq_threshold.
  # top_n: if 'all', plot all motifs found. If top_n is an integer, then all top top_n motifs are plotted.
  # plot_curves: if TRUE, plot all the curves with coloured motifs.
  
  ### check input ############################################################################################
  # check freq_threshold
  oldpar <- par(no.readonly = TRUE)
  # Ensure the original settings are restored when the function exits
  on.exit(par(oldpar))
  if(max(motifs_search_results$V_frequencies)<freq_threshold)
    stop('There are no motifs with frequency at least equal to freq_threshold.')
  # check top_n
  if(top_n!='all'){
    if(length(top_n)!=1)
      stop('top_n not valid.')
    if(top_n%%1!=0)
      stop('top_n should be an integer.')
    if(top_n<1)
      stop('top_n should be at least 1.')
  }
  
  ### select motifs to plot ##################################################################################
  d=ncol(motifs_search_results$Y0[[1]])
  N=length(motifs_search_results$Y0)
  index_plot=which(motifs_search_results$V_frequencies>=freq_threshold)
  if(top_n!='all'){
    if(length(index_plot)>top_n)
      index_plot=index_plot[seq_len(top_n)]
  }
  K=length(index_plot)
  V0=motifs_search_results$V0[index_plot]
  V1=motifs_search_results$V1[index_plot]
  V_dom=lapply(V0,function(v) rowSums(!is.na(v))!=0)
  V_length=motifs_search_results$V_length[index_plot]
  V_occurrences=motifs_search_results$V_occurrences[index_plot]
  V_frequencies=motifs_search_results$V_frequencies[index_plot]
  V_mean_diss=motifs_search_results$V_mean_diss[index_plot]
  R_motifs=motifs_search_results$R_motifs[index_plot]
  ### plot motifs ############################################################################################
  layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
  lapply(seq_len(d),
         function(j){
           par(mar=c(3,4,4,2)+0.1)
           plot(V0[[1]][,j],type='l',col=rainbow(K),lwd=5,lty=1,main=paste0(ylab,"-","Dimension:",j),xlim=c(1,max(V_length)),
                ylab=ylab,ylim=c(min(unlist(V0)),max(unlist(V0))))
           mapply(function(v,k) points(v[,j],type='l',col=rainbow(K)[k+1],lwd=5,lty=1,ylab=ylab),
                  V0[-1],seq_len(K-1))
           par(mar=c(0,0,0,0))
           return()})
  plot.new()
  legend('left',paste('motif',seq_len(K)),col=rainbow(K),lwd=7,lty=1,bty="n",xpd=TRUE)
  if(!is.null(V1[[1]])){
    layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             plot(V1[[1]][,j],type='l',col=rainbow(K),lwd=5,lty=1,main=paste(ylab,"-","Dimension:",j,'derivative'),xlim=c(1,max(V_length)),
                  ylab=ylab,ylim=c(min(unlist(V1)),max(unlist(V1))))
             mapply(function(v,k) points(v[,j],type='l',col=rainbow(K)[k+1],lwd=5,lty=1,ylab=ylab),
                    V1[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             return()})
    plot.new()
    legend('left',paste('motif',seq_len(K)),col=rainbow(K),lwd=7,lty=1,bty="n",xpd=TRUE)
  }
  
  ### plot motifs with matched curves ########################################################################
  if(is.null(V1[[1]])){
    mapply(function(v,v_dom,v_occurrences,v_frequencies,k,R_motif){
      Y_inters_k=mapply(function(y,s_k_i,v_dom){
        v_len=length(v_dom)
        Y_inters_k=as.matrix(as.matrix(y[s_k_i-1+seq_len(v_len),])[v_dom,])
        return(Y_inters_k)},  # c('curve','shift','diss')
        motifs_search_results$Y0[v_occurrences[,1]],v_occurrences[,2],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
      layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
      Y0_diff_k=lapply(Y_inters_k,
                       function(Y0_inters_k){
                         y0_min=apply(Y0_inters_k, 2, min, na.rm = TRUE)
                         y0_max=apply(Y0_inters_k, 2, max, na.rm = TRUE)
                         y0_diff=y0_max-y0_min
                         return(y0_diff)
                       })
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
               if(transformed){
                 y_plot[v_dom,]=Reduce('cbind',
                                       mapply(function(Y_inters_k, Y_diff_k) {
                                         y0_min=min(Y_inters_k[,j])
                                         y0_norm = t( (t(Y_inters_k[,j]) - y0_min) / Y_diff_k[j] )
                                         y0_const = (Y_diff_k[j] == 0)
                                         y0_norm[,y0_const] = 0.5
                                         return(y0_norm)},
                                         Y_inters_k, Y0_diff_k, SIMPLIFY=FALSE) )
               } else {
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
               }
               matplot(y_plot,type='l',col=v_occurrences[,1]+1,lwd=round(-4/R_motif*v_occurrences[,3]+5,2),
                       lty=1,ylab=ylab,main=paste0('Motif ',k,' (',v_frequencies,' occurrences) - ',ylab,"-","Dimension:",j))
               points(v[,j],type='l',col='black',lwd=7,lty=1)
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
             })
      return()},V0,V_dom,V_occurrences,V_frequencies,seq_len(K),R_motifs)
  }else{
    mapply(function(v0,v1,v_dom,v_occurrences,v_frequencies,k,R_motif){
      Y0_inters_k=mapply(
        function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          Y_inters_k=as.matrix(as.matrix(y[s_k_i-1+seq_len(v_len),])[v_dom,])
          return(Y_inters_k)},
        motifs_search_results$Y0[v_occurrences[,1]],v_occurrences[,2],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
      Y1_inters_k=mapply(
        function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          Y_inters_k=as.matrix(as.matrix(y[s_k_i-1+seq_len(v_len),])[v_dom,])
          return(Y_inters_k)},
        motifs_search_results$Y1[v_occurrences[,1]],v_occurrences[,2],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
      layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(7,1))
      Y0_diff_k=lapply(Y0_inters_k,
                       function(Y0_inters_k){
                         y0_min=apply(Y0_inters_k, 2, min, na.rm = TRUE)
                         y0_max=apply(Y0_inters_k, 2, max, na.rm = TRUE)
                         y0_diff=y0_max-y0_min
                         return(y0_diff)
                       })
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y0_inters_k))
               if(transformed){
                 y_plot[v_dom,]=Reduce('cbind',
                                       mapply(function(Y_inters_k, Y_diff_k) {
                                         y0_min=min(Y_inters_k[,j])
                                         y0_norm = t( (t(Y_inters_k[,j]) - y0_min[j]) / Y_diff_k[j] )
                                         y0_const = (Y_diff_k[j] == 0)
                                         y0_norm[,y0_const] = 0.5
                                         return(y0_norm)},
                                         Y0_inters_k, Y0_diff_k, SIMPLIFY=FALSE) 
                 )
               } else {
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))
               }
               matplot(y_plot,type='l',col=v_occurrences[,1]+1,lwd=round(-4/R_motif*v_occurrences[,3]+5,2),
                       lty=1,ylab=ylab,main=paste0('Motif ',k,' (',v_frequencies,' occurrences) - ',ylab,"-","Dimension:",j))
               points(v0[,j],type='l',col='black',lwd=7,lty=1)
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
             })
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y1_inters_k))
               if(transformed){
                 y_plot[v_dom,]=Reduce('cbind',
                                       mapply(function(Y1_inters_k, Y_diff_k) {
                                         y1_norm = t( t(Y1_inters_k[,j])/ Y_diff_k[j] )
                                         y0_const = (Y_diff_k[j] == 0)
                                         y1_norm[,y0_const] = 0
                                         return(y1_norm)},
                                         Y1_inters_k, Y0_diff_k, SIMPLIFY=FALSE)
                 ) 
               } else {
                 y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))
               }               
               matplot(y_plot,type='l',col=v_occurrences[,1]+1,lwd=round(-4/R_motif*v_occurrences[,3]+5,2),
                       lty=1,ylab=ylab,main=paste0('Motif ',k,' (',v_frequencies,' occurrences) - ',ylab,"-","Dimension:",j,' derivative'))
               points(v1[,j],type='l',col='black',lwd=7,lty=1)
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
             })
      return()},V0,V1,V_dom,V_occurrences,V_frequencies,seq_len(K),R_motifs)
  }
  
  ### plot curves with motifs ################################################################################
  if(plot_curves){
    if(is.null(motifs_search_results$Y1[[1]])){
      mapply(function(y0,i){
        s_i=lapply(V_occurrences,function(occurrences) occurrences[occurrences[,1]==i,2])
        motifs_in_curve=rep(seq_len(K),unlist(lapply(s_i,length)))
        s_i=unlist(s_i)
        Y_inters_k=mapply(function(v_dom,s_i_k,y){
          v_len=length(v_dom)
          Y_inters_k=matrix(NA,nrow=length(v_dom),ncol=d)
          Y_inters_k[v_dom,]=as.matrix(as.matrix(y[s_i_k-1+seq_len(v_len),])[v_dom,])
          return(Y_inters_k)},
          V_dom[motifs_in_curve],s_i,MoreArgs=list(y0),SIMPLIFY=FALSE)
        layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 plot(y0[,j],type='l',main=paste('Region',i,'-',ylab,"-","Dimension:",j),ylab=ylab)
                 for(k in seq_along(motifs_in_curve)){
                   lines(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]),Y_inters_k[[k]][,j],col=rainbow(K)[motifs_in_curve[k]],lwd=5)
                   rect(s_i[k], min(y0[,j],na.rm=TRUE)-10, tail(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]), n=1), max(y0[,j],na.rm=TRUE)+10,
                        border = scales::alpha(rainbow(K)[motifs_in_curve[k]], 0.05), col = scales::alpha(rainbow(K)[motifs_in_curve[k]], 0.05))
                 }
               })
        plot.new()
        if(length(motifs_in_curve)==0){
          legend_text=''
        }else{
          legend_text=paste('motif',unique(motifs_in_curve))
        }
        legend('left',legend_text,col=rainbow(K)[unique(motifs_in_curve)],lwd=7,lty=1,bty="n",xpd=TRUE,title='Motifs')
        return()},motifs_search_results$Y0,seq_len(N))
    }else{
      mapply(function(y0,y1,i){
        s_i=lapply(V_occurrences,function(occurrences) occurrences[occurrences[,1]==i,2])
        motifs_in_curve=rep(seq_len(K),unlist(lapply(s_i,length)))
        s_i=unlist(s_i)
        Y0_inters_k=mapply(function(v_dom,s_i_k,y){
          v_len=length(v_dom)
          Y_inters_k=matrix(NA,nrow=length(v_dom),ncol=d)
          Y_inters_k[v_dom,]=as.matrix(as.matrix(y[s_i_k-1+seq_len(v_len),])[v_dom,])
          return(Y_inters_k)},
          V_dom[motifs_in_curve],s_i,MoreArgs=list(y0),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(function(v_dom,s_i_k,y){
          v_len=length(v_dom)
          Y_inters_k=matrix(NA,nrow=length(v_dom),ncol=d)
          Y_inters_k[v_dom,]=as.matrix(as.matrix(y[s_i_k-1+seq_len(v_len),])[v_dom,])
          return(Y_inters_k)},
          V_dom[motifs_in_curve],s_i,MoreArgs=list(y1),SIMPLIFY=FALSE)
        layout(matrix(c(seq_len(d),rep(d+1,d)),ncol=2),widths=c(7,1))
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 plot(y0[,j],type='l',main=paste('Region',i,'-',ylab,"-","Dimension:",j),ylab=ylab,xlab='')
                 for(k in seq_along(motifs_in_curve)){
                   lines(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]),Y0_inters_k[[k]][,j],col=rainbow(K)[motifs_in_curve[k]],lwd=5)
                 }
               })
        plot.new()
        if(length(motifs_in_curve)==0){
          legend_text=''
        }else{
          legend_text=paste('motif',unique(motifs_in_curve))
        }
        legend('left',legend_text,col=rainbow(K)[unique(motifs_in_curve)],lwd=7,lty=1,bty="n",xpd=TRUE,title='Motifs')
        lapply(seq_len(d),
               function(j){
                 par(mar=c(3,4,4,2)+0.1)
                 plot(y1[,j],type='l',main=paste('Region',i,'-',paste(ylab,"-","Dimension:",j,'derivative')),ylab=ylab)
                 for(k in seq_along(motifs_in_curve)){
                   lines(s_i[k]-1+seq_len(V_length[motifs_in_curve[k]]),Y1_inters_k[[k]][,j],col=rainbow(K)[motifs_in_curve[k]],lwd=5)
                 }
               })
        plot.new()
        if(length(motifs_in_curve)==0){
          legend_text=''
        }else{
          legend_text=paste('motif',unique(motifs_in_curve))
        }
        legend('left',legend_text,col=rainbow(K)[unique(motifs_in_curve)],lwd=7,lty=1,bty="n",xpd=TRUE,title='Motifs')
        return()},motifs_search_results$Y0,motifs_search_results$Y1,seq_len(N))
    }
  }
  
  return()
}