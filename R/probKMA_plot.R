#' @title Plot the Results of probKMA
#'
#' @description
#' The `probKMA_plot` function visualizes the results obtained from the `probKMA` analysis. It generates a series of plots
#' including motif memberships across different curves, the progression of the objective function over iterations,
#' and the Bhattacharyya distance between memberships. Depending on the parameters, it can plot both original
#' and cleaned motifs across multiple dimensions, providing insights into the embedding and characteristics of
#' identified motifs.
#'
#' @param probKMA_results A list containing the output from the `probKMA` function. This list should include elements such as:
#'   \describe{
#'     \item{Y0}{A list of matrices representing the original curves.}
#'     \item{Y1}{A list of matrices representing the derivatives of the curves (if applicable).}
#'     \item{V0}{A list of motifs.}
#'     \item{V1}{A list of derived motifs (if applicable).}
#'     \item{P}{A matrix indicating motif memberships across curves.}
#'     \item{S}{A matrix indicating motif start positions in curves.}
#'     \item{S_clean}{A matrix indicating cleaned motif start positions (if applicable).}
#'     \item{P_clean}{A matrix indicating cleaned motif memberships (if applicable).}
#'     \item{V0_clean}{A list of cleaned motifs (if applicable).}
#'     \item{V1_clean}{A list of derived cleaned motifs (if applicable).}
#'     \item{iter}{An integer indicating the number of iterations performed.}
#'     \item{J_iter}{A numeric vector recording the objective function value at each iteration.}
#'     \item{BC_dist_iter}{A numeric vector recording the Bhattacharyya distance between memberships at each iteration.}
#'   }
#'
#' @param plot A logical flag indicating whether to produce the plots. If `TRUE`, the function generates all relevant plots. If `FALSE`, no plots are produced.
#'
#' @param ylab A character vector of length `d`, providing labels for the y-axis in each dimension. Defaults to an empty string (`''`) for all dimensions.
#'
#' @param sil_avg A numeric vector containing the average silhouette scores for each embedded motif. This parameter is used to annotate the plots with silhouette information. Defaults to `NULL`, meaning no silhouette scores are displayed.
#'
#' @param cleaned A logical value indicating whether to plot only the cleaned motifs (`TRUE`) or all motifs (`FALSE`). When set to `TRUE`, the function highlights motifs that have been cleaned based on predefined criteria. Defaults to `FALSE`.
#'
#' @param transformed A logical value indicating whether to normalize the curve segments to the interval [0,1] before applying the dissimilarity measure. Setting `transformed = TRUE` scales each curve segment between 0 and 1, which allows for the identification of motifs with consistent shapes but different amplitudes. This normalization is useful for cases where motif occurrences may vary in amplitude but have similar shapes, enabling better pattern recognition across diverse data scales.
#'
#' @return 
#' The function generates a series of plots visualizing:
#' \itemize{
#'   \item \strong{Motifs with Matched Curves}: Displays the original curves with embedded motifs overlaid. If `cleaned = TRUE`, only cleaned motifs are highlighted.
#'   \item \strong{Memberships}: Shows bar plots representing the membership scores of each motif across all curves.
#'   \item \strong{Objective Function}: Plots the progression of the objective function (`J_iter`) over iterations to illustrate convergence.
#'   \item \strong{Bhattacharyya Distance}: Plots the Bhattacharyya distance (`BC_dist_iter`) between memberships over iterations to assess similarity.
#' }
#' No value is returned; the function is used solely for its side effects of generating visualizations.
#'
#' @details
#' The `probKMA_plot` function performs the following operations:
#' \enumerate{
#'   \item **Motif Visualization**:
#'     \itemize{
#'       \item Plots the original curves (`Y0`) with embedded motifs (`V0`). If derivatives (`Y1` and `V1`) are available, additional plots are generated for them.
#'       \item When `cleaned = TRUE`, the function highlights only the cleaned motifs (`V0_clean` and `V1_clean`), providing a clearer view of significant motifs.
#'       \item Utilizes color coding and legends to differentiate between different motifs and their instances.
#'     }
#'   \item **Memberships**:
#'     \itemize{
#'       \item Generates bar plots showing the membership scores (`P` or `P_clean`) of each motif across all curves.
#'       \item Provides a visual representation of how strongly each motif is associated with different curves.
#'     }
#'   \item **Objective Function and Bhattacharyya Distance**:
#'     \itemize{
#'       \item Plots the objective function (`J_iter`) over the iterations to demonstrate the optimization process.
#'       \item Plots the Bhattacharyya distance (`BC_dist_iter`) to measure the similarity between motif memberships across iterations.
#'     }
#' }
#' 
#' The function is designed to handle multiple dimensions (`d`) and can accommodate both original and derivative data if provided. 
#' It also supports the visualization of cleaned motifs, which are motifs that have been refined based on specific criteria to ensure quality and relevance.
#'
#' @export
probKMA_plot <- function(probKMA_results, plot, ylab='',sil_avg=NULL,cleaned = FALSE, transformed = FALSE){
  d=ncol(probKMA_results$Y0[[1]])
  N=nrow(probKMA_results$P)
  K=ncol(probKMA_results$P)
  V_dom=lapply(probKMA_results$V0,function(v) rowSums(!is.na(v))!=0)
  S_k=split(probKMA_results$S,rep(seq_len(K),each=N))
  P_k=split(probKMA_results$P,rep(seq_len(K),each=N))
  
  ### plot motifs with matched curves #######################################################################
  if(cleaned){
    S_clean_k=split(probKMA_results$S_clean,rep(seq_len(K),each=N))
    P_clean_k=split(probKMA_results$P_clean,rep(seq_len(K),each=N))
    S_clean_i=split(probKMA_results$S_clean,seq_len(nrow(probKMA_results$S_clean)))
    P_clean_i=split(probKMA_results$P_clean,seq_len(nrow(probKMA_results$S_clean)))
    has_a_motif <- sapply(P_clean_i,function(p_i){sum(p_i)!=0}) # at least one motif embedded
    if(is.null(probKMA_results$V1[[1]])){
      NULL
      #### Plot curves all together with the embedded motifs -------------
      # colorful plot (deleted)
      # mapply(function(v0,v_dom,s_k,p_clean_k,k)
      #   {
      #   layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
      #   keep=which(p_clean_k==1)
      #   # matrix extension (eventually)
      #   full_mat <- probKMA_results$Y0
      #   # max len of the curves
      #   maxLen <-  max(sapply(full_mat, nrow))
      #   # Extended full matrix
      #   full_mat <- sapply(full_mat, padding, maxLen,simplify = FALSE)
      #   # Join
      #   full_mat <- do.call(cbind, full_mat)
      #   dom_sequence <- sapply(s_k[keep], function(x){seq(x,x+max(0,length(v_dom)-1))},simplify=TRUE)
      #   # Extract keep submatrix
      #   column <- numeric(length(keep) * d)
      #   for (i in 1:length(keep)) {
      #     start_col <- (keep[i] - 1) * d + 1
      #     end_col <- start_col + d - 1
      #     column[((i - 1) * d + 1):((i - 1) * d + d)] <- start_col:end_col
      #   }
      #   # Plot curves all together 
      #   lapply(seq_len(d),
      #          function(j){
      #            par(mar=c(3,4,4,2)+0.1)
      #            matplot(full_mat[,seq(j,ncol(full_mat),by = d)],col=scales::alpha('gray30',0.30),lwd=1.5,lty=1,type = 'l',
      #                    main=paste0('Curves & ','Motif_',k,' - Dimension:',j,'\n',
      #                                'Number of instances: ',ncol(dom_sequence),' - sil_avg: ',
      #                                ifelse(is.null(sil_avg),"",round(sil_avg[k], digits = 3))))
      #            univariate_mat <- as.matrix(as.matrix(full_mat[,column])[,seq(j,d*length(keep),by = d)])
      #            matplot(univariate_mat,type="l",col=rainbow(length(keep)),lwd=1.5,lty=5,
      #                    main=paste0('Curves & ','Motif_',k,' - Dimension:',j,'\n',
      #                                'Number of instances: ',ncol(dom_sequence),' - sil_avg: ',
      #                                ifelse(is.null(sil_avg),"",round(sil_avg[k], digits = 3))),add=TRUE)
      #            for (i in 1:ncol(dom_sequence)) {
      #              lines(dom_sequence[,i],univariate_mat[dom_sequence[,i],i], col = 'red',lwd=2.5)
      #              rect(dom_sequence[1,i], min(univariate_mat,na.rm=TRUE)-10, tail(dom_sequence[,i], n=1), max(univariate_mat,na.rm=TRUE)+10,
      #                   border = scales::alpha(rainbow(length(keep))[i], 0.05), col = scales::alpha(rainbow(length(keep))[i], 0.05))
      #            }
      #            title(sub = paste('Curves & ', 'Motif:', k, ' - Dimension:', j))
      #            par(mar=c(0,0,0,0))
      #            plot.new()
      #            legend("topright", legend = c(paste0("c",(1:N)[keep]),paste0('motif_',k)),
      #                   col = c(rainbow(length(keep)),'red'), lty = 1,lwd=c(rep(1,length(keep)),4),cex=1.5)
      #          })
      #   return()},probKMA_results$V0_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
      # 
      #### for each curve plot the embedded motifs  ############################################################################################
      # single line plot
      
      # mapply(function(curve_i,s_i,p_clean_i,i){
      #   layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
      #   lapply(seq_len(d),
      #          function(j)
      #            {
      #             par(mar=c(3,4,4,2)+0.1)
      #             univariate_mat <- curve_i[,j]
      #             matplot(univariate_mat,type="l",col='black',lwd=1.0,lty=1,
      #                     main=paste0('c',i,' - Dimension:',j),
      #                     #ylab=expression(partialdiff * paste(ylab, "/", partialdiff, "t")),
      #                     xlab='domain')
      #             # plot motifs with different colors 
      #             col_i <- 1
      #             col <- rainbow(sum(p_clean_i))
      #             lapply(seq_len(K),function(k)
      #              {
      #               if(p_clean_i[k]==1)
      #               {
      #                 univariate_motif <- probKMA_results$V0_clean[[k]][,j]
      #                 v_dom <- V_dom[[k]]
      #                 dom_seq <- seq.int(s_i[k],s_i[k]+max(0,length(v_dom)-1))
      #                 if(transformed)
      #                 {
      #                   max_num = max(univariate_mat[dom_seq[v_dom]], na.rm=TRUE)
      #                   min_num = min(univariate_mat[dom_seq[v_dom]], na.rm=TRUE)
      #                   univariate_motif[v_dom] = univariate_motif[v_dom]*(max_num - min_num) + min_num
      #                 }
      #                 lines(dom_seq[v_dom], univariate_motif[v_dom], col = col[col_i] ,lwd=2.5)
      #                 rect(dom_seq[v_dom][1],min(univariate_mat,na.rm=TRUE)-10, tail(dom_seq[v_dom], n=1), max(univariate_mat,na.rm=TRUE)+10,
      #                                 border = scales::alpha(col[col_i], 0.1), col = scales::alpha(col[col_i], 0.1)) 
      #                 col_i <<- col_i + 1
      #               }
      #               })
      #               par(mar=c(0,0,0,0))
      #               plot.new()
      #               legend('left',legend=paste0('motif',which(p_clean_i==1)),col=col,lwd=7,lty=1,bty="n",xpd=TRUE)
      #               })
      #   return()},probKMA_results$Y0[ has_a_motif] # consider curves with at least one embedded motif
      #             ,S_clean_i[ has_a_motif],P_clean_i[has_a_motif],seq_len(N)[has_a_motif])
      # mapply(function(v,v_dom,s_k,p_clean_k,k){
      #   layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
      #   keep=which(p_clean_k==1)
      #   Y_inters_k=mapply(
      #     function(y,s_k_i,v_dom){
      #       v_len=length(v_dom)
      #       d=ncol(y)
      #       y_len=nrow(y)
      #       index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
      #       Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
      #                        matrix(y[index[index<=y_len],],ncol=d),
      #                        matrix(NA,nrow=sum(index>y_len),ncol=d))
      #       return(Y_inters_k)},
      #     probKMA_results$Y0[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
      #   Y0_diff_k=lapply(Y_inters_k,
      #                    function(Y0_inters_k){
      #                      y0_min=apply(Y0_inters_k, 2, min, na.rm = TRUE)
      #                      y0_max=apply(Y0_inters_k, 2, max, na.rm = TRUE)
      #                      y0_diff=y0_max-y0_min
      #                      return(y0_diff)
      #                    })
      #   lapply(seq_len(d),
      #          function(j){
      #            par(mar=c(3,4,4,2)+0.1)
      #            y_plot=matrix(NA,nrow=length(v_dom),ncol=length(Y_inters_k))
      #            if(transformed){
      #              y_plot[v_dom,]=Reduce('cbind',
      #                                    mapply(function(Y_inters_k, Y_diff_k) {
      #                                      y0_min=min(Y_inters_k[,j])
      #                                      y0_norm = t( (t(Y_inters_k[,j]) - y0_min) / Y_diff_k[j] )
      #                                      y0_const = (Y_diff_k[j] == 0)
      #                                      y0_norm[,y0_const] = 0.5
      #                                      return(y0_norm)},
      #                                      Y_inters_k, Y0_diff_k, SIMPLIFY=FALSE) ) 
        #          }else{
        #            y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))
        #          }
        #          
        #          matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
        #          points(v[,j],type='l',col='black',lwd=7,lty=1)
        #          par(mar=c(0,0,0,0))
        #          plot.new()
        #          legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
        #        })
        # return()},
        # probKMA_results$V0_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
    }else{
      print("working on it")
      ##### Plot curves and derivatives all together with the embedded motifs  ############################################################################################
      # colorful plot with everything together
      
      # mapply(function(v0,v1,v_dom,s_k,p_clean_k,k)
      # {
      #   layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
      #   keep=which(p_clean_k==1)
      #   # matrix extension (eventually)
      #   full_mat <- probKMA_results$Y0
      #   full_mat_dev <- probKMA_results$Y1
      #   # max len of the curves
      #   maxLen <-  max(sapply(full_mat, nrow))
      #   # Extended full matrix
      #   full_mat <- sapply(full_mat, padding, maxLen,simplify = FALSE)
      #   full_mat_dev <- sapply(full_mat_dev, padding, maxLen,simplify = FALSE)
      #   # Join
      #   full_mat <- do.call(cbind, full_mat)
      #   full_mat_dev <- do.call(cbind, full_mat_dev)
      #   dom_sequence <- sapply(s_k[keep], function(x){seq(x,x+max(0,length(v_dom)-1))},simplify=TRUE)
      #   # Extract "keep" submatrix
      #   column <- numeric(length(keep) * d)
      #   for (i in 1:length(keep)) {
      #     start_col <- (keep[i] - 1) * d + 1
      #     end_col <- start_col + d - 1
      #     column[((i - 1) * d + 1):((i - 1) * d + d)] <- start_col:end_col
      #   }
      #   # Plot curves all together 
      #   lapply(seq_len(d),
      #          function(j){
      #            par(mar=c(3,4,4,2)+0.1)
      #            matplot(full_mat[,seq(j,ncol(full_mat),by = d)],col=scales::alpha('gray30',0.30),lwd=1.5,lty=1,type = 'l',
      #                    main=paste0('Curves & ','Motif_',k,' - Dimension:',j,'\n',
      #                                'Number of instances: ',ncol(dom_sequence),' - sil_avg: ',
      #                                ifelse(is.null(sil_avg),"",round(sil_avg[k], digits = 3))))
      #            univariate_mat <- as.matrix(as.matrix(full_mat[,column])[,seq(j,d*length(keep),by = d)])
      #            matplot(univariate_mat,type="l",col=rainbow(length(keep)),lwd=1.5,lty=5,
      #                    main=paste0('Curves & ','Motif_',k,' - Dimension:',j,'\n',
      #                                'Number of instances: ',ncol(dom_sequence),' - sil_avg: ',
      #                                ifelse(is.null(sil_avg),"",round(sil_avg[k], digits = 3))),add = TRUE)
      #            for (i in 1:ncol(dom_sequence)) {
      #              lines(dom_sequence[,i],univariate_mat[dom_sequence[,i],i], col = 'red',lwd=2.5)
      #              rect(dom_sequence[1,i], min(univariate_mat,na.rm=TRUE)-10, tail(dom_sequence[,i], n=1), max(univariate_mat,na.rm=TRUE)+10,
      #                   border = scales::alpha(rainbow(length(keep))[i], 0.05), col = scales::alpha(rainbow(length(keep))[i], 0.05))
      #            }
      #            title(sub = paste('Curves & ', 'Motif_', k, ' - Dimension:', j))
      #            par(mar=c(0,0,0,0))
      #            plot.new()
      #            legend("topright", legend = c(paste0("c",(1:N)[keep]),paste0('motif_',k)),
      #                   col = c(rainbow(length(keep)),'red'), lty = 1,lwd=c(rep(1,length(keep)),4),cex=1.5)
      #          })
      #   # Plot Derivatives all together 
      #   lapply(seq_len(d),
      #          function(j){
      #            par(mar=c(3,4,4,2)+0.1)
      #            matplot(full_mat_dev[,seq(j,ncol(full_mat_dev),by = d)],col=scales::alpha('gray30',0.30),lwd=1.5,lty=1,type = 'l',
      #                    main=paste0('Derivatives & ','Motif_',k,' - Dimension:',j,'\n',
      #                                'Number of instances: ',ncol(dom_sequence),' - sil_avg: ',
      #                                ifelse(is.null(sil_avg),"",round(sil_avg[k], digits = 3))))
      #            univariate_mat <- as.matrix(as.matrix(full_mat_dev[,column])[,seq(j,length(keep),by = d)])
      #            matplot(univariate_mat,type="l",col=rainbow(length(keep)),lwd=1.5,lty=5,
      #                    main=paste0('Derivatives & ','Motif_',k,' - Dimension:',j,'\n',
      #                                'Number of instances: ',ncol(dom_sequence),' - sil_avg: ',
      #                                ifelse(is.null(sil_avg),"",round(sil_avg[k], digits = 3))),add = TRUE)
      #            for (i in 1:ncol(dom_sequence)) {
      #              lines(dom_sequence[,i],univariate_mat[dom_sequence[,i],i], col = 'red',lwd=2.5)
      #              rect(dom_sequence[1,i], min(univariate_mat,na.rm=TRUE)-10, tail(dom_sequence[,i], n=1), max(univariate_mat,na.rm=TRUE)+10,
      #                   border = scales::alpha(rainbow(length(keep))[i], 0.05), col = scales::alpha(rainbow(length(keep))[i], 0.05))
      #            }
      #            title(sub = paste('Derivatives & ', 'Motif:', k, ' - Dimension:', j))
      #            par(mar=c(0,0,0,0))
      #            plot.new()
      #            legend("topright", legend = c(paste0("c",(1:N)[keep]),paste0('motif_',k)),
      #                   col = c(rainbow(length(keep)),'red'), lty = 1,lwd=c(rep(1,length(keep)),4),cex=1.5)
      #          })
      #   
      #   return()},probKMA_results$V0_clean,probKMA_results$V1_clean,
      #             V_dom,S_clean_k,P_clean_k,seq_len(K))
      
    ##### for each curve plot the embedded motifs  ############################################################################################
    # for each curve the motif highlighted
      
      # mapply(function(curve_i,dev_i,s_i,p_clean_i,i){
      # #plot curve with embedded motif 
      # layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
      # lapply(seq_len(d),
      #        function(j)
      #        {
      #          par(mar=c(3,4,4,2)+0.1)
      #          univariate_mat <- curve_i[,j]
      #          matplot(univariate_mat,type="l",col='black',lwd=1.0,lty=1,
      #                  main=paste0('c',i,' - Dimension:',j),
      #                  #ylab=expression(partialdiff * paste(ylab, "/", partialdiff, "t")),
      #                  xlab='domain')
      #          # plot motifs with different colors 
      #          col_i <- 1
      #          col <- rainbow(sum(p_clean_i))
      #          lapply(seq_len(K),function(k)
      #          {
      #            if(p_clean_i[k]==1)
      #            {
      #              univariate_motif <- probKMA_results$V0_clean[[k]][,j]
      #              v_dom <- V_dom[[k]]
      #              dom_seq <- seq.int(s_i[k],s_i[k]+max(0,length(v_dom)-1))
      #              if(transformed)
      #              {
      #                max_num = max(univariate_mat[dom_seq[v_dom]], na.rm=TRUE)
      #                min_num = min(univariate_mat[dom_seq[v_dom]], na.rm=TRUE)
      #                univariate_motif[v_dom] = univariate_motif[v_dom]*(max_num - min_num) + min_num
      #              }
      #              lines(dom_seq[v_dom], univariate_motif[v_dom], col = col[col_i] ,lwd=2.5)
      #              rect(dom_seq[v_dom][1],min(univariate_mat,na.rm=TRUE)-10, tail(dom_seq[v_dom], n=1), max(univariate_mat,na.rm=TRUE)+10,
      #                   border = scales::alpha(col[col_i], 0.1), col = scales::alpha(col[col_i], 0.1)) 
      #              col_i <<- col_i + 1
      #            }
      #          })
      #          par(mar=c(0,0,0,0))
      #          plot.new()
      #          legend('left',legend=paste0('motif_',which(p_clean_i==1)),col=col,lwd=7,lty=1,bty="n",xpd=TRUE)
      #        })
      
      #plot derivative with embedded motif 
      #for each derivative curve the motif highlighted
      
      # layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
      # lapply(seq_len(d),
      #        function(j)
      #        {
      #          par(mar=c(3,4,4,2)+0.1)
      #          univariate_mat_curve <- curve_i[,j]
      #          univariate_mat <- dev_i[,j]
      #          matplot(univariate_mat,type="l",col='black',lwd=1.0,lty=1,
      #                  main=paste0('c',i,' - Dimension:',j),
      #                  ylab=expression(partialdiff * paste(ylab, "/", partialdiff, "t")),
      #                  xlab='domain')
      #          # plot motifs with different colors 
      #          col_i <- 1
      #          col <- rainbow(sum(p_clean_i))
      #          lapply(seq_len(K),function(k)
      #          {
      #            if(p_clean_i[k]==1)
      #            {
      #              univariate_motif <- probKMA_results$V1_clean[[k]][,j]
      #              v_dom <- V_dom[[k]]
      #              dom_seq <- seq.int(s_i[k],s_i[k]+max(0,length(v_dom)-1))
      #              if(transformed)
      #              {
      #                max_num = max(univariate_mat_curve[dom_seq[v_dom]], na.rm=TRUE)
      #                min_num = min(univariate_mat_curve[dom_seq[v_dom]], na.rm=TRUE)
      #                univariate_motif[v_dom] = univariate_motif[v_dom]*(max_num - min_num)
      #              }
      #              lines(dom_seq[v_dom], univariate_motif[v_dom], col = col[col_i] ,lwd=2.5)
      #              rect(dom_seq[v_dom][1],min(univariate_mat,na.rm=TRUE)-10, tail(dom_seq[v_dom], n=1), max(univariate_mat,na.rm=TRUE)+10,
      #                   border = scales::alpha(col[col_i], 0.1), col = scales::alpha(col[col_i], 0.1)) 
      #              col_i <<- col_i + 1
      #            }
      #          })
      #          par(mar=c(0,0,0,0))
      #          plot.new()
      #          legend('left',legend=paste0('motif_',which(p_clean_i==1)),col=col,lwd=7,lty=1,bty="n",xpd=TRUE)
      #        })
      # 
      # return()},probKMA_results$Y0[ has_a_motif], # consider curves with at least one embedded motif
      #           probKMA_results$Y1[ has_a_motif],
      #           S_clean_i[ has_a_motif],P_clean_i[has_a_motif],seq_len(N)[has_a_motif])
    
      mapply(function(v0,v1,v_dom,s_k,p_clean_k,k){
        keep=which(p_clean_k==1)
        Y0_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          probKMA_results$Y0[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(
          function(y,s_k_i,v_dom){
            v_len=length(v_dom)
            d=ncol(y)
            y_len=nrow(y)
            index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
            Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                             matrix(y[index[index<=y_len],],ncol=d),
                             matrix(NA,nrow=sum(index>y_len),ncol=d))
            return(Y_inters_k)},
          probKMA_results$Y1[keep],s_k[keep],MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
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
                                           y0_norm = t( (t(Y_inters_k[,j]) - y0_min) / Y_diff_k[j] )
                                           y0_const = (Y_diff_k[j] == 0)
                                           y0_norm[,y0_const] = 0.5
                                           return(y0_norm)},
                                           Y0_inters_k, Y0_diff_k, SIMPLIFY=FALSE) )
                 }else{
                   y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))}
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
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
                 }else{
                   y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))}
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=1,lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j],' derivative'))
                 points(v1[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0_clean,probKMA_results$V1_clean,V_dom,S_clean_k,P_clean_k,seq_len(K))
    }
    
  }
  else
  {
    if(is.null(probKMA_results$V1[[1]])){
      mapply(function(v,v_dom,s_k,p_k,k){
        Y_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          probKMA_results$Y0,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
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
                 }else{
                   y_plot[v_dom,]=Reduce('cbind',lapply(Y_inters_k,function(Y_inters_k) Y_inters_k[,j]))}
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
                 points(v[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0,V_dom,S_k,P_k,seq_len(K))
    }else{
      mapply(function(v0,v1,v_dom,s_k,p_k,k){
        Y0_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          probKMA_results$Y0,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        Y1_inters_k=mapply(function(y,s_k_i,v_dom){
          v_len=length(v_dom)
          d=ncol(y)
          y_len=nrow(y)
          index=max(1,s_k_i)-1+seq_len(v_len-max(0,1-s_k_i))
          Y_inters_k=rbind(matrix(NA,nrow=max(0,1-s_k_i),ncol=d),
                           matrix(y[index[index<=y_len],],ncol=d),
                           matrix(NA,nrow=sum(index>y_len),ncol=d))
          return(Y_inters_k)},
          probKMA_results$Y1,s_k,MoreArgs=list(v_dom),SIMPLIFY=FALSE)
        layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
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
                                           y0_norm = t( (t(Y_inters_k[,j]) - y0_min) / Y_diff_k[j] )
                                           y0_const = (Y_diff_k[j] == 0)
                                           y0_norm[,y0_const] = 0.5
                                           return(y0_norm)},
                                           Y0_inters_k, Y0_diff_k, SIMPLIFY=FALSE) )
                 }else{
                   y_plot[v_dom,]=Reduce('cbind',lapply(Y0_inters_k,function(Y_inters_k) Y_inters_k[,j]))}
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j]))
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
                 }else{
                   y_plot[v_dom,]=Reduce('cbind',lapply(Y1_inters_k,function(Y_inters_k) Y_inters_k[,j]))}
                 matplot(y_plot,type='l',col=seq_len(N)+1,lwd=round(5*p_k,2),lty=1,ylab=ylab[j],main=paste('Motif',k,'-',ylab[j],' derivative'))
                 points(v1[,j],type='l',col='black',lwd=7,lty=1)
                 par(mar=c(0,0,0,0))
                 plot.new()
                 legend('left',legend='motif center',col='black',lwd=7,lty=1,bty="n",xpd=TRUE)
               })
        return()},
        probKMA_results$V0,probKMA_results$V1,V_dom,S_k,P_k,seq_len(K))
    }
  }
  
  layout(matrix(1:(2*d), ncol=2, byrow=TRUE), widths=c(8.5,1))
  ### plot motifs ############################################################################################
  if(cleaned){
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             motif_length=unlist(lapply(probKMA_results$V0_clean,nrow))
             plot(probKMA_results$V0_clean[[1]][,j],type='l',col=2,lwd=7,lty=1,main=ylab[j],xlim=c(1,max(motif_length)),
                  ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V0_clean),na.rm=TRUE),max(unlist(probKMA_results$V0_clean),na.rm=TRUE)))
             mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                    probKMA_results$V0_clean[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             plot.new()
             legend('left',paste('motif_',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
             return()})
  }else{
    lapply(seq_len(d),
           function(j){
             par(mar=c(3,4,4,2)+0.1)
             motif_length=unlist(lapply(probKMA_results$V0,nrow))
             plot(probKMA_results$V0[[1]][,j],type='l',col=2,lwd=7,lty=1,main=ylab[j],xlim=c(1,max(motif_length)),
                  ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V0),na.rm=TRUE),max(unlist(probKMA_results$V0),na.rm=TRUE)))
             mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                    probKMA_results$V0[-1],seq_len(K-1))
             par(mar=c(0,0,0,0))
             plot.new()
             legend('left',paste('motif_',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
             return()})
  }
  if(!is.null(probKMA_results$V1[[1]])){
    layout(matrix(1:(2*d),ncol=2,byrow=TRUE),widths=c(8.5,1))
    if(cleaned){
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               motif_length=unlist(lapply(probKMA_results$V1_clean,nrow))
               plot(probKMA_results$V1_clean[[1]][,j],type='l',col=2,lwd=7,lty=1,main=paste0(ylab[j],' derivative'),xlim=c(1,max(motif_length)),
                    ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V1_clean),na.rm=TRUE),max(unlist(probKMA_results$V1_clean),na.rm=TRUE)))
               mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                      probKMA_results$V1_clean[-1],seq_len(K-1))
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',paste('motif_',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
               return()})
    }else{
      lapply(seq_len(d),
             function(j){
               par(mar=c(3,4,4,2)+0.1)
               motif_length=unlist(lapply(probKMA_results$V1,nrow))
               plot(probKMA_results$V1[[1]][,j],type='l',col=2,lwd=7,lty=1,main=paste0(ylab[j],' derivative'),xlim=c(1,max(motif_length)),
                    ylab=ylab[j],ylim=c(min(unlist(probKMA_results$V1),na.rm=TRUE),max(unlist(probKMA_results$V1),na.rm=TRUE)))
               mapply(function(v,k) points(v[,j],type='l',col=k+2,lwd=7,lty=1,ylab=ylab),
                      probKMA_results$V1[-1],seq_len(K-1))
               par(mar=c(0,0,0,0))
               plot.new()
               legend('left',paste('motif_',seq_len(K)),col=seq_len(K)+1,lwd=7,lty=1,bty="n",xpd=TRUE)
               return()})
    }
  }
  
  ### plot memberships #######################################################################################
  par(mfrow=c(K,1),mar=c(3,4,4,2)+0.1)
  if(cleaned){
    mapply(function(p_k,p_clean_k,k){
      col=rep('lightgray',N)
      col[p_clean_k==1]='gray35'
      barplot(p_k,names.arg=seq_len(N),col=col,las=2,ylim=c(0,1),ylab='memberships',main=paste('Motif_',k))
      legend('left',paste('p_clean==1'),col='gray35',pch=15,cex=2,bty="n",xpd=TRUE)
    },P_k,P_clean_k,seq_len(K))
  }else{
    mapply(function(p_k,k){
      barplot(p_k,names.arg=seq_len(N),col='gray',las=2,ylim=c(0,1),ylab='memberships',main=paste('Motif_',k))
    },P_k,seq_len(K))
  }
  
  ### plot objective function and Bhattacharyya distance #####################################################
  par(mfrow=c(1,1))
  plot(seq_len(probKMA_results$iter),probKMA_results$J_iter,type='l',xlab='iter',ylab='objective function',
       main=paste0('Objective function Jm:',round(tail(probKMA_results$J_iter,n=1)),4),lwd=2)
  points(seq_len(probKMA_results$iter), probKMA_results$J_iter,
         col = 'red', pch = 8, cex = 2,lwd=3)
  for (i in 1:probKMA_results$iter) {
    segments(i,0.0,i,probKMA_results$J_iter[i], col = "red",lwd=0.65)
  }
  
  plot(seq_len(probKMA_results$iter),probKMA_results$BC_dist_iter,type='l',xlab='iter',ylab='distance between memberships',
       main=paste0('Bhattacharyya distance between memberships',round(tail(probKMA_results$BC_dist_iter,n=1)),4),lwd=2)
  points(seq_len(probKMA_results$iter), probKMA_results$BC_dist_iter,
         col = 'red', pch = 8, cex = 2,lwd=3)
  for (i in 1:probKMA_results$iter) {
    segments(i,0.0,i,probKMA_results$BC_dist_iter[i], col = "red",lwd=0.65)
  }
  return()
}  
