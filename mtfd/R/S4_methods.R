#' @title showDetails
#' @description Shows details of an object of class MyS4Class.
#' @param object An object of class MyS4Class.
#' @export
# error_str can be a vector or a matrix with n_cols = motif_len
setGeneric("generateCurves", function(object,error_type,error_str) {
  standardGeneric("generateCurves")
})

#' @rdname showDetails
#' @export
setMethod("generateCurves", "motifSimulation", function(object,error_type,error_str) {
  breaks=lapply(rep(object@len,length.out=object@N),seq,from=0,by=object@dist_knots) # generate N lists with equally spaced nodes 'dist_knots' from 0 to 'len'
  basis=lapply(breaks,function(breaks_i) create.bspline.basis(norder=object@norder,breaks=breaks_i)) # nbasis = norder + nbreaks - 2 = norder + interior_knots,
  len_motifs <- unlist(lapply(object@mot_details,function(x){x$len}))
  # loop for each curve
  fd_curves <- NULL
  if(error_type == 'coeff') {
    tryCatch({
      # Attempt to coerce the input to a numeric vector
      error_str <- as.numeric(error_str)
    }, error = function(e) {
      stop("Input cannot be vectorized.")
    })
    fd_curves <- mapply(function(motifs_in_curves_i,len_i,basis_i,j,len_motifs){
      list_coeff <- NULL
      if(is.numeric(object@distribution)) {
        or_coeff=sample(object@distribution,size = len_i,replace = TRUE)
      } else if(object@distribution=='unif'){
        or_coeff=runif(len_i,min=object@coeff_min,max=object@coeff_max) # coefficients of the curve = degree of freedom are initialized uniformly
      } else if(object@distribution=='beta') {  
        or_coeff=object@coeff_min+rbeta(len_i,0.45,0.45)*(object@coeff_max-object@coeff_min) # coefficients of the curve = degree of freedom are initialized uniformly
      } else {
        stop('Wrong \'distrib\'')
      }
      coeff <- or_coeff
      fda_no_error=fda::fd(coef=coeff,basisobj=basis_i)
      or_y_no_error <- mtfd:::generate_curve_vector(fda_no_error)
      if(!is.null(motifs_in_curves_i)) {
        pos_coeff_motifs=unlist(mapply(
          function(a,b) seq(a)+b,
          rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
          motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
        list_coeff <- lapply(error_str,function(sd_noise) {
          coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
            print(paste(" --- Adding motif", id, "to curve", j,"with noise ",sd_noise))
            object@mot_details[[id]]$weights})) +
            rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
          # For each chosen coefficient add a uniform number and a gaussian noise
          return(coeff)
        })
      }
      fda_with_error <- NULL
      or_y <- NULL
      
      if(!is.null(motifs_in_curves_i)) {
        pos_coeff_motifs=unlist(mapply(
          function(a,b) seq(a)+b,
          rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
          motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE))
        coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
                                                 object@mot_details[[id]]$weights}))
        fda_with_motif <- fda::fd(coef=coeff,basisobj=basis_i)
        or_y_motif <- mtfd:::generate_curve_vector(fda_with_motif)
        
        fda_with_error= Map(fda::fd,list_coeff,MoreArgs = list(basisobj=basis_i)) # Fitting curves using such coefficients and basis
        or_y <- Map(mtfd:::generate_curve_vector,fda_with_error)
        num_motifs <- length(motifs_in_curves_i$motif_id)
        SNR <- c()
        for (k in seq_along(error_str)) {
          SNR_num <- numeric(num_motifs)
          SNR_den <- numeric(num_motifs)
          
          # Loop over each motif
          for (n in seq_along(motifs_in_curves_i$motif_id)) {
            start_break <- motifs_in_curves_i$starting_coeff_pos[n]
            start_point <- (start_break - 1) * object@dist_knots
            end_point   <- start_point + len_motifs[n]
            
            # Calculate variance for SNR numerator and denominator
            SNR_num[n] <- var(or_y_no_error[start_point:end_point])
            SNR_den[n] <- var(or_y[[k]][start_point:end_point] - or_y_no_error[start_point:end_point])
          }
          # Calculate SNR in decibels for the k-th error
          SNR <- c(SNR,10 * log10(SNR_num / SNR_den))
        }
        return(list(basis = basis_i,
                    background = list(or_coeff = or_coeff,no_error_y = or_y_no_error),
                    no_noise = list(or_coeff = coeff,motif_y = or_y_motif),
                    with_error = list(error_structure = error_str,error_y = or_y),SNR = SNR))
      }
      return(list(basis = basis_i,
                  background = list(or_coeff = or_coeff,no_error_y = or_y_no_error),
                  with_error = list(error_structure = error_str,error_y = or_y)))
    },object@motifs_in_curves,object@norder-1+rep(object@len/object@dist_knots,length.out=object@N),basis,1:object@N,MoreArgs = list(len_motifs),SIMPLIFY=FALSE)
    
  }else if(error_type == 'pointwise'){
    fd_curves <- lapply(1:object@N, function(x){
      coeff <- NULL
      len_i <- object@norder-1+object@len/object@dist_knots
      if(is.numeric(object@distribution)) {
        coeff=sample(object@distribution,size = len_i,replace = TRUE)
      }
      else if(object@distribution=='unif'){
        coeff=runif(len_i,min=object@coeff_min,max=object@coeff_max) # If we don't have motifs we sample len_i coefficients uniformely
      }else if(object@distribution=='beta'){
        coeff=object@coeff_min+rbeta(len_i,0.45,0.45)*(object@coeff_max-object@coeff_min) 
      }else {
        stop('Wrong \'distrib\'')
      }
      mtfd:::generate_background_curve(object@len, object@dist_knots, object@norder,coeff, add_noise = TRUE)
    })
    for(i in 1:length(object@mot_details)){
      print(paste("Dealing with motif", i))
      temp_mot <- object@mot_details[[i]]
      curve_ids <- unique(temp_mot$appearance$curve)
      for(j in curve_ids){
        print(paste(" --- Adding motif", i, "to curve", j))
        temp_curve <- fd_curves[[j]]
        temp_pattern <- (temp_mot$appearance %>% filter(curve == j)) %>% dplyr::select(motif_id, start_break_pos)
        fd_curves[[j]] <- mtfd:::add_motif(base_curve  = temp_curve,
                                    mot_pattern = temp_pattern,
                                    mot_len     = temp_mot$len,
                                    dist_knots  = object@dist_knots,
                                    mot_order   = object@norder,
                                    mot_weights = temp_mot$weight,
                                    error_str   = error_str)
      }
    }
    fd_curves <- mtfd:::.transform_list(fd_curves,error_str)
  } else {
    stop("\'error_type\' must be choosen between \'coeff\' and \'pointwise\'")
  }
  return(fd_curves = fd_curves)
})

#' @export
setGeneric("plot", function(object,curves,path) 
  standardGeneric("plot")
)

#' @export
setMethod("plot",c(object = "motifSimulation", curves = "list", path = "character"),function(object,curves,path) {
  output_file <- file.path(path, "plots.pdf")
  
  # Create the directory if it does not exist
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Open a PDF device with the correct path
  pdf(file = output_file, width = 8, height = 6) 
  for (k in seq_along(curves)) {
    if(!is.null(curves[[k]]$with_error$error_y)) {
      curve_data_no_error <- data.frame(
        t = seq(0, curves[[k]]$basis$rangeval[2]-1),
        x = curves[[k]]$background$no_error_y,
        z = curves[[k]]$no_noise$motif_y)
      names(curve_data_no_error) <- c("t","x","z")
 
      curve_data_error <- NULL
      curve_data_error <- data.frame(
        t = seq(0, curves[[k]]$basis$rangeval[2]-1),
        x = curves[[k]]$with_error$error_y,
        SNR = unlist(rep(curves[[k]]$SNR[1:length(curves[[k]]$SNR)],len = curves[[k]]$basis$rangeval[2])))
      names(curve_data_error) <- c("t",paste0("x",seq(length(curves[[k]]$with_error$error_y))),"SNR")
      
      motif_lines <- mapply(function(id_motif, pos_motif, instance) {
        motif_t = seq((pos_motif - 1) * object@dist_knots,
                      (pos_motif - 1) * object@dist_knots + object@mot_details[[id_motif]]$len)
        motif_x = lapply(curves[[k]]$with_error$error_y,function(curve){return(curve[motif_t + 1])})
        
        return(lapply(motif_x,function(motif){ data.frame(t = motif_t, x = motif, motif_id = factor(paste(id_motif, instance, sep = "_")),
                                                          initial_number = str_extract(as.character(id_motif), "^[^_]+"),
                                                          xmin = (pos_motif - 1) * object@dist_knots,
                                                          xmax = (pos_motif - 1) * object@dist_knots + object@mot_details[[id_motif]]$len)}))
      }, object@motifs_in_curves[[k]]$motif_id, object@motifs_in_curves[[k]]$starting_coeff_pos, seq_along(object@motifs_in_curves[[k]]$motif_id), SIMPLIFY = FALSE)
      
      motif_colors <- c( "1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange",
                         "5" = "purple", "6" = "cyan", "7" = "magenta", "8" = "brown",
                         "9" = "pink", "10" = "grey")
      motif_colors <- rep(motif_colors,length.out = length(object@mot_details))
      if(length(object@mot_details) > 10 )
        attr(motif_colors,"names")[11:length(object@mot_details)] <- as.character(as.integer(attr(motif_colors,"names")[11:length(object@mot_details)]) + 10)
      
      max_dataframes <- max(sapply(motif_lines, function(sublist) length(sublist)))
      # Inizializza una lista per memorizzare i risultati
      motif_data <- vector("list", max_dataframes)
      # Cicla su ogni "livello" dei data frame
      for (i in seq_len(max_dataframes)) {
        # Estrai i data frame dal livello i-esimo
        dataframes_at_level_i <- lapply(motif_lines, function(sublist) {
          sublist[[i]]
        })
        # Fai il bind_rows sui data frame estratti
        motif_data[[i]] <- bind_rows(dataframes_at_level_i)
        names(motif_data[[i]]) <- c("t","x","motif_id","initial_number","xmin","xmax")
      }
      p <- lapply(1:length(motif_data), function(j) {
        # Creare i label dei motif_id
        motif_labels <- paste("motif_id:", unique(motif_data[[j]]$initial_number))
        
        pic <- ggplot() +
          # Plot the main curve in gray30 and add it to the legend
          geom_line(data = curve_data_no_error, aes(x = t, y = x, color = 'background_curve'), linewidth = 0.5) +
          # Plot the curve with motifs in gold and add it to the legend
          geom_line(data = curve_data_no_error, aes(x = t, y = z, color = 'with_motif'), linewidth = 0.5) +
          # Plot the error curve in black
          geom_line(data = curve_data_error, aes_string(x = "t", y = paste0("x", j)), color = "black", linewidth = 0.5) +
          # Add shaded rectangles for motif positions with transparency
          geom_rect(data = motif_data[[j]], 
                    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = factor(initial_number)), 
                    alpha = 0.005) +
          # Add SNR text on top of each rectangle
          geom_text(data = motif_data[[j]], 
                    aes(x = (xmin + xmax) / 2, y = Inf, 
                        label = paste("SNR:", round(curve_data_error$SNR[j], 3))), 
                    vjust = 1.5, color = "black", size = 3.5) +  
          # Plot motifs with distinct colors
          geom_line(data = motif_data[[j]], aes(x = t, y = x, color = factor(initial_number), group = motif_id), linewidth = 1.0) + 
          # Add color and fill scales with custom labels for motif_id
          scale_color_manual(
            values = c('background_curve' = scales::alpha('gray30', 0.15), 
                       'with_motif' = 'gold', 
                       motif_colors),
            labels = c('background_curve' = 'background_curve', 
                       'with_motif' = 'with_motif', 
                       setNames(motif_labels, unique(motif_data[[j]]$initial_number)))  # Correctly map labels to motif ids
          ) +
          scale_fill_manual(
            values = motif_colors, 
            labels = setNames(motif_labels, unique(motif_data[[j]]$initial_number))  # Correctly map labels to motif ids
          ) +
          # Title without SNR
          labs(
            title = paste0('<b><span style="color:#0073C2;">Random curve ', k,' - type_error ',ifelse(j%%length(curves[[k]]$with_error$error_y)==0,length(curves[[k]]$with_error$error_y),j%%length(curves[[k]]$with_error$error_y)),'</span></b>'), 
            x = "t", 
            y = paste0("x", j)
          ) +
          # Apply a clean theme with a subtle grid
          theme_minimal(base_size = 15) +
          theme(
            plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),  # Enhanced title with colors
            axis.title = element_text(size = 14, margin = margin(t = 10)),  # Add margin to axis titles
            axis.text = element_text(size = 12),  # Reduce axis text size
            legend.text = element_text(size = 12),  # Adjust legend text size
            legend.position = "right",  # Position the legend to the right
            legend.box.margin = margin(10, 10, 10, 10),  # Add margin around the legend
            plot.margin = margin(15, 15, 15, 15),  # Increase plot margins
            panel.grid.major = element_line(color = "gray90"),  # Lighten the grid lines
            panel.grid.minor = element_blank()  # Remove minor grid lines for a cleaner look
          ) +
          guides(
            color = guide_legend(ncol = 1, byrow = TRUE, title = NULL),  # Remove the title from the color legend
            fill = FALSE  # Hide the fill legend
          )
        
        return(pic)
      })
      
      Map(print,p)
    }
  }
  dev.off()
})
