#' @title showDetails
#' @description Shows details of an object of class MyS4Class.
#' @param object An object of class MyS4Class.
#' @export
# noise_str can be a vector or a matrix with n_cols = motif_len
setGeneric("generateCurves", function(object,noise_type = NULL, noise_str = NULL,seed_background = 777,seed_motif = 43213,
                                      only_der = TRUE,coeff_min_shift = -10,coeff_max_shift = 10) {
  standardGeneric("generateCurves")
})

#' @rdname showDetails
#' @export
setMethod("generateCurves", "motifSimulation", function(object,noise_type = NULL, noise_str = NULL,seed_background = 777,seed_motif = 43213,
                                                        only_der = TRUE,coeff_min_shift = -10,coeff_max_shift = 10) {
  breaks=lapply(rep(object@len,length.out=object@N),seq,from=0,by=object@dist_knots) # generate N lists with equally spaced nodes 'dist_knots' from 0 to 'len'
  basis=lapply(breaks,function(breaks_i) create.bspline.basis(norder=object@norder,breaks=breaks_i)) # nbasis = norder + nbreaks - 2 = norder + interior_knots,
  len_motifs <- unlist(lapply(object@mot_details,function(x){x$len}))
  set.seed(seed_background)
  # loop for each curve
  fd_curves <- NULL
  if(is.null(noise_type) && !is.null(noise_str)) {
    stop("\'noise_str'\ must be specified because \'noise_type\' is not null")
  }
  if(is.null(noise_str) && !is.null(noise_type)) {
    stop("\'noise_type'\ must be specified because \'noise_str\' is not null")
  }
  if((is.null(noise_type) && !is_empty(object@motifs_in_curves))) {
    stop("\'noise_type\' must be choosen between \'coeff\' and \'pointiwse\'")
  } 
  if(is.null(noise_type)) {
    noise_type <- "pointwise" # to generate background curves
  }
  if(!is.null(noise_type) && noise_type == 'coeff') {
    tryCatch({
      # Attempt to coerce the input to a numeric vector
      noise_str <- lapply(noise_str,function(error){return(as.numeric(error))})
    }, error = function(e) {
      stop("\'noise_str'\ cannot be vectorized.")
    })
    if(length(noise_str) != length(object@mot_details)) {
     stop("\'noise_str\' must have the same length of \'mot_details\'") 
    }
    fd_curves <- mapply(function(motifs_in_curves_i,len_i,basis_i,j,len_motifs){
      list_coeff_error <- NULL
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
      fda_with_error <- NULL
      or_y <- NULL
      shifted_coeff <- NULL
      set.seed(seed_motif)
      if(!is.null(motifs_in_curves_i)) {
        pos_coeff_motifs=unlist(mapply(
          function(a,b) seq(a)+b,
          rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
          motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
        shifted_coeff <- rep(runif(length(motifs_in_curves_i$motif_id),min=coeff_min_shift,max=coeff_max_shift),rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1) 
        for(z in 1:length(noise_str[[1]])) {
        list_coeff_error <- lapply(1:length(noise_str[[1]]),function(null) {
          sd_noise <- NULL
          if(only_der) {
          coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
            sd_noise <<- noise_str[[id]][z]
            print(paste(" --- Adding motif", id, "to curve", j,"with noise ",sd_noise))
            object@mot_details[[id]]$coeffs})) +
            rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
          } else {
            coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
              sd_noise <<- noise_str[[id]][z]
              print(paste(" --- Adding motif", id, "to curve", j,"with noise ",sd_noise))
              object@mot_details[[id]]$coeffs})) + shifted_coeff + 
                rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
          }
          # For each chosen coefficient add a uniform number and a gaussian noise
          return(coeff)
        })
        }
      if(only_der) {
        coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
                                  object@mot_details[[id]]$coeffs}))
      } else {
        coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
                                    object@mot_details[[id]]$coeffs})) + shifted_coeff
      }
      fda_with_motif <- fda::fd(coef=coeff,basisobj=basis_i)
      or_y_motif <- mtfd:::generate_curve_vector(fda_with_motif)
      
      fda_with_error <- Map(fda::fd,list_coeff_error,MoreArgs = list(basisobj=basis_i)) # Fitting curves using such coefficients and basis
      or_y <- Map(mtfd:::generate_curve_vector,fda_with_error)
      num_motifs <- length(motifs_in_curves_i$motif_id)
      SNR <- vector("list",length(noise_str[[1]]))
      for (k in seq_along(noise_str[[1]])) {
        SNR_num <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
        SNR_den <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
        
        # Loop over each motif
        for (n in seq_along(motifs_in_curves_i$motif_id)) {
          start_break <- motifs_in_curves_i$starting_coeff_pos[n]
          start_point <- (start_break - 1) * object@dist_knots
          end_point   <- start_point + len_motifs[n]
          
          # Calculate variance for SNR numerator and denominator
          SNR_num <- rbind(SNR_num,data.frame(xmin =start_point,xmax = end_point,SNR = var(or_y_motif[start_point:end_point])))
          SNR_den <- rbind(SNR_den,data.frame(xmin =start_point,xmax = end_point,SNR = var(or_y[[k]][start_point:end_point] - or_y_motif[start_point:end_point])))
        }
        # Calculate SNR in decibels for the k-th error
        SNR[[k]] <- SNR_num
        SNR[[k]]$SNR <- 10 * log10(SNR_num$SNR / SNR_den$SNR) #transform SNR in decibel
      }
      return(list(basis = basis_i,
                  background = list(or_coeff = or_coeff,no_error_y = or_y_no_error),
                  no_noise = list(or_coeff = coeff,motif_y = or_y_motif),
                  with_noise = list(noise_structure = noise_str,noise_y = or_y),
                  SNR = SNR))
      }
        
      return(list(basis = basis_i,
                  background = list(or_coeff = or_coeff,no_error_y = or_y_no_error)))
    },object@motifs_in_curves,object@norder-1+rep(object@len/object@dist_knots,length.out=object@N),basis,1:object@N,MoreArgs = list(len_motifs),SIMPLIFY=FALSE)
    
  }else if(!is.null(noise_type) && noise_type == 'pointwise'){
    if(length(noise_str) != length(object@mot_details)) {
      stop("\'noise_str\' must have the same length of \'mot_details\'") 
    }
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
    set.seed(seed_motif)
    curve_ids <- unique(unlist(sapply(object@mot_details,function(mot_details){mot_details$occurrences$curve})))
    for(j in curve_ids){
      print(paste(" --- Adding motifs to curve", j))
      temp_curve <- fd_curves[[j]]
      temp_pattern <- do.call(rbind,lapply(object@mot_details,function(mot_details){(mot_details$occurrences %>% filter(curve == j)) %>% dplyr::select(motif_id, start_break_pos)}))
      temp_len <- do.call(rbind,lapply(unique(apply(temp_pattern,1,function(row){row[1]})),function(k){data.frame("motif_id" = k,"len" = object@mot_details[[k]]$len)}))
      temp_weights <- lapply(as.character(unique(apply(temp_pattern, 1, function(row) { row[1] }))), 
                             function(k) { 
                               object@mot_details[[as.numeric(k)]]$coeffs
                             })
      names(temp_weights) <- as.character(unique(apply(temp_pattern, 1, function(row) { row[1] })))
      fd_curves[[j]] <- mtfd:::add_motif(base_curve  = temp_curve,
                                  mot_pattern = temp_pattern,
                                  mot_len     = temp_len,
                                  dist_knots  = object@dist_knots,
                                  mot_order   = object@norder,
                                  mot_weights = temp_weights,
                                  noise_str   = noise_str,
                                  only_der = only_der,
                                  coeff_min_shift = coeff_min_shift,
                                  coeff_max_shift = coeff_max_shift)
    }
    fd_curves <- mtfd:::.transform_list(fd_curves,noise_str)
  } else if(!is.null(noise_type)){
    stop("\'noise_type\' must be choosen between \'coeff\' and \'pointwise\'")
  }
  return(fd_curves = fd_curves)
})

#' @export
setGeneric("plot", function(object,curves,path) 
  standardGeneric("plot")
)

#' @export
setMethod("plot",c(object = "motifSimulation", curves = "list", path = "character"),
          function(object,curves,path) {
  output_file <- file.path(path, "plots.pdf")
  
  # Create the directory if it does not exist
  if (!dir.exists(path)) {
    dir.create(path)
  }
  # Open a PDF device with the correct path
  pdf(file = output_file, width = 8, height = 6) 
  for (k in seq_along(curves)) {
    if(purrr::is_empty(object@motifs_in_curves)) {
      curve_data_no_error <- data.frame(
        t = seq(0, curves[[k]]$basis$rangeval[2]-1),
        x = curves[[k]]$background$no_error_y
      )
      names(curve_data_no_error) <- c("t", "x")
      
      p <- ggplot() +
        # Plot the main curve in gray30
        geom_line(data = curve_data_no_error, aes(x = t, y = x, color = 'background_curve'), linewidth = 0.5) +
        
        scale_color_manual(
          values = c('background_curve' = scales::alpha('gray30', 0.90))
        ) +
        
        # Title
        labs(
          title = paste0('<b><span style="color:#0073C2;">Random curve ', k, '</span></b>'),
          x = "t",
          y = "x"
        ) +
        
        # Clean theme with subtle grid
        theme_minimal(base_size = 15) +
        theme(
          plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
          axis.title = element_text(size = 14, margin = margin(t = 10)),
          axis.text = element_text(size = 12),
          legend.position = "none", 
          plot.margin = margin(15, 15, 15, 15),  # Margini intorno al grafico
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank()
        ) 
      
      print(p)
    }else if(!is.null(curves[[k]]$with_noise)) {
      curve_data_no_error <- data.frame(
        t = seq(0, curves[[k]]$basis$rangeval[2]-1),
        x = curves[[k]]$background$no_error_y,
        z = curves[[k]]$no_noise$motif_y)
      names(curve_data_no_error) <- c("t","x","z")
 
      curve_data_error <- NULL
      curve_data_error <- data.frame(
        t = seq(0, curves[[k]]$basis$rangeval[2]-1),
        x = curves[[k]]$with_noise$noise_y)
      names(curve_data_error) <- c("t",paste0("x",seq(length(curves[[k]]$with_noise$noise_y))))
      
      motif_lines <- mapply(function(id_motif, pos_motif, instance) {
        motif_t = seq((pos_motif - 1) * object@dist_knots,
                      (pos_motif - 1) * object@dist_knots + object@mot_details[[id_motif]]$len)
        motif_x = lapply(curves[[k]]$with_noise$noise_y,function(curve){return(curve[motif_t + 1])})
        
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
      # Initialize a list to store results
      motif_data <- vector("list", max_dataframes)
      for (i in seq_len(max_dataframes)) {
        # Extract data frames at the i-th level
        dataframes_at_level_i <- lapply(motif_lines, function(sublist) sublist[[i]])
        # Combine the extracted data frames
        motif_data[[i]] <- bind_rows(dataframes_at_level_i)
        names(motif_data[[i]]) <- c("t", "x", "motif_id", "initial_number", "xmin", "xmax")
      }
      p <- lapply(1:length(motif_data), function(j) {
        # Create motif_id labels
        motif_labels <- paste("motif_id:", unique(motif_data[[j]]$initial_number))
        pic <- ggplot() +
          # Plot the main curve in gray30
          geom_line(data = curve_data_no_error, aes(x = t, y = x, color = 'background_curve'), linewidth = 0.5) +
          # Plot the curve with motifs in gold
          geom_line(data = curve_data_no_error, aes(x = t, y = z, color = 'zero_noise_motif'), linewidth = 0.5) +
          # Plot the error curve
          geom_line(data = curve_data_error, aes_string(x = "t", y = paste0("x", j)), color = "black", linewidth = 0.5) +
          # Add shaded rectangles for motif positions
          geom_rect(data = motif_data[[j]], 
                    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = factor(initial_number)), 
                    alpha = 0.005) +
          # Add SNR text on top of each rectangle
          geom_text(data = curves[[k]]$SNR[[j]], 
                    aes(x = (xmin + xmax) / 2, y = Inf, 
                        label = paste("SNR:", round(SNR, 3))),
                    vjust = 1.5, color = "black", size = 3.5) +
          # Plot motifs with distinct colors
          geom_line(data = motif_data[[j]], aes(x = t, y = x, color = factor(initial_number), group = motif_id), linewidth = 1.0) + 
          # Add color and fill scales with custom labels for motif_id
          scale_color_manual(
            values = c('background_curve' = scales::alpha('gray30', 0.15), 
                       'zero_noise_motif' = 'gold', 
                       motif_colors),
            labels = c('background_curve' = 'background_curve', 
                       'zero_noise_motif' = 'zero_noise_motif', 
                       setNames(motif_labels, unique(motif_data[[j]]$initial_number)))
          ) +
          scale_fill_manual(
            values = motif_colors, 
            labels = setNames(motif_labels, unique(motif_data[[j]]$initial_number))
          ) +
          # Title
          labs(
            title = paste0('<b><span style="color:#0073C2;">Random curve ', k, ' - type_error ', 
                           ifelse(j %% length(curves[[k]]$with_noise$noise_y) == 0, 
                                  length(curves[[k]]$with_noise$noise_y), j %% length(curves[[k]]$with_noise$noise_y)), 
                           '</span></b>'), 
            x = "t", 
            y = "x"
          ) +
          # Clean theme with subtle grid
          theme_minimal(base_size = 15) +
          theme(
            plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
            axis.title = element_text(size = 14, margin = margin(t = 10)),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.position = "right",
            legend.box.margin = margin(10, 10, 10, 10),
            plot.margin = margin(15, 15, 15, 15),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank()
          ) +
          guides(
            color = guide_legend(ncol = 1, byrow = TRUE, title = NULL),
            fill = "none",
          )
        
        return(pic)
      })
      
      Map(print, p)
    }
  }
  if(!purrr::is_empty(object@motifs_in_curves)) {
    curves_data <- list()
    for(motif_id in 1:length(object@mot_details)) {
      for(curve_k in 1:length(object@motifs_in_curves)) {
        motif_instance <- 1
        if(!is.null(object@motifs_in_curves[[curve_k]])) {
          for(z in 1:length(object@motifs_in_curves[[curve_k]]$motif_id)) {
            if(object@motifs_in_curves[[curve_k]]$motif_id[z] == motif_id) {
              # Calcola l'intervallo della curva
              start <- (object@motifs_in_curves[[curve_k]]$starting_coeff_pos[z] - 1) * object@dist_knots + 1
              
              # Crea un identificatore unico per ogni istanza
              instance_id <- paste0("curve_", curve_k, "_", motif_instance)
              motif_instance <- motif_instance + 1
              
              curve_y <- curves[[curve_k]]$no_noise$motif_y[start:(start + object@mot_details[[motif_id]]$len - 1)]
              
              # Creazione di una sequenza x per il plotting
              x <- seq_along(curve_y)
              
              # Aggiungi i dati alla lista
              curves_data[[length(curves_data) + 1]] <- data.frame(x = x, y = curve_y, id = motif_id, instance = instance_id)  
            }
          }
        }
      } 
    }
  
    # Combina i dati in un unico dataframe
    curves_df <- bind_rows(curves_data)
    
    # Ottieni la lista unica degli ID
    unique_ids <- unique(curves_df$id)
    # Loop per plottare ogni ID separatamente su pagine diverse
    for (id in unique_ids) {
      # Filtra i dati per il singolo ID
      plot_data <- curves_df %>% filter(id == !!id)
      
      # Crea l'oggetto fd e il dataframe per motif_y
      mot_breaks <- seq(from = 0, to = length(curve_y), by = object@dist_knots)
      mot_basis  <- create.bspline.basis(norder = object@norder, breaks = mot_breaks)
      fd_curve <- fd(coef = object@mot_details[[id]]$coeffs, basisobj = mot_basis)
      motif_y <- mtfd:::generate_curve_vector(fd_curve = fd_curve)
      motif_y <- data.frame(x = seq_along(motif_y), y = motif_y, instance = "motif_y")
      names(motif_y) <- c("x","y","instance")
      p <- ggplot() +
        geom_line(data = plot_data, aes(x = x, y = y, color = instance), size = 1, linetype = "longdash") +
        
        geom_line(data = motif_y, aes(x = x, y = y, color = instance), size = 2, linetype = "solid") +  
        
        scale_color_manual(
          values = c("motif_y" = "black",setNames(rainbow(length(unique(plot_data$instance))), unique(plot_data$instance))),
          labels = c("motif_y" = paste0("motif ", id),unique(plot_data$instance))
        ) +
        labs(
          title = paste0('<b><span style="color:#0073C2;">Motif ', id, '</span></b>'), 
          x = "X-axis", 
          y = "Y-axis",
          color = "Legend Title"  
        ) +
        theme_minimal(base_size = 15) +
        theme(
          plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
          axis.title = element_text(size = 14, margin = margin(t = 10)),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "right",
          legend.box.margin = margin(10, 10, 10, 10),
          plot.margin = margin(15, 15, 15, 15),
          panel.grid.major = element_line(color = "gray90"),
          panel.grid.minor = element_blank()
        ) +
        guides(
          color = guide_legend(ncol = 1, byrow = TRUE, title = NULL),
          fill = "none"
        )
        
        print(p)
    }
  }
  dev.off()
})


#' @export
setGeneric("to_motifDiscovery", function(curves) 
  standardGeneric("to_motifDiscovery")
)

#' @export
setMethod("to_motifDiscovery",c(curves = "list"), function(curves) {
  
  if(purrr::is_empty(curves)) {
    warning("\'curves' object is empty. Return \'NULL\'")
    return(NULL)
  }
  
  n <- 1
  for(i in 1:length(curves)) {
    if("with_noise" %in% names(curves[[i]])) {
      n <- length(curves[[i]]$with_noise$noise_y)
      break
    }
  }
 
  result <- vector("list", n)  
  for(i in 1:n) {
    result[[i]] <- setNames(lapply(curves, function(curve) {  
      if("with_noise" %in% names(curve)) {
        return(curve$with_noise$noise_y[[i]])
      } else {
        return(curve$background$no_error_y)
      }
    }), paste0("c", seq_along(curves)))  
  }
  
  result <- setNames(result, paste0("Noise_type_", seq_len(n)))
  
  return(result)
})
