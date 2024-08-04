#' @title showDetails
#' @description Shows details of an object of class MyS4Class.
#' @param object An object of class MyS4Class.
#' @export
# error_str can be a vector or a matrix with n_cols = motif_len
setGeneric("generateCurves", function(object,error_str) {
  standardGeneric("generateCurves")
})

#' @rdname showDetails
#' @export
setMethod("generateCurves", "motifSimulation", function(object,error_str) {
  breaks=lapply(rep(object@len,length.out=object@N),seq,from=0,by=object@dist_knots) # generate N lists with equally spaced nodes 'dist_knots' from 0 to 'len'
  basis=lapply(breaks,function(breaks_i) create.bspline.basis(norder=object@norder,breaks=breaks_i)) # nbasis = norder + nbreaks - 2 = norder + interior_knots,
  len_motifs <- unlist(lapply(object@mot_details,function(x){x$len}))
  # loop for each curve
  fd_curves <- NULL
  if(!object@is_appearance_defined) {
    if(is.data.frame(error_str) || is.matrix(error_str)) {
      stop("error_str must be a single number or a vector")
    }
    fd_curves <- mapply(function(motifs_in_curves_i,len_i,basis_i,j,len_motifs){
      list_coeff <- NULL
      if(object@distribution=='unif'){
        or_coeff=runif(len_i,min=object@coeff_min,max=object@coeff_max) # coefficients of the curve = degree of freedom are initialized uniformly
        coeff <- or_coeff
        fda_no_error=fda::fd(coef=coeff,basisobj=basis_i)
        or_y_no_error <- mtfd:::generate_curve_vector(fda_no_error)
        if(!is.null(motifs_in_curves_i)) {
          pos_coeff_motifs=unlist(mapply(
            function(a,b) seq(a)+b,
            rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
            motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
          list_coeff <- lapply(error_str,function(sd_noise) {
            coeff[pos_coeff_motifs]= unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
              print(paste(" --- Adding motif", id, "to curve", j,"with noise ",sd_noise))
              object@mot_details[[id]]$weights})) +
              rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
            # For each chosen coefficient add a uniform number and a gaussian noise
            coeff[-pos_coeff_motifs]=runif(len_i-length(pos_coeff_motifs),min=object@coeff_min,max=object@coeff_max) # All other curve coefficients are randomly generated
            return(coeff)
          })
        }
        # Same as before but sampling from a beta
      }else if(object@distribution=='beta'){
        or_coeff=object@coeff_min+rbeta(len_i,0.45,0.45)*(object@coeff_max-object@coeff_min) # coefficients of the curve = degree of freedom are initialized uniformly
        coeff <- or_coeff
        fda_no_error=fda::fd(coef=coeff,basisobj=basis_i)
        or_y_no_error <- mtfd:::generate_curve_vector(fda_no_error)
        if(!is.null(motifs_in_curves_i)) {
          pos_coeff_motifs=unlist(mapply(
            function(a,b) seq(a)+b,
            rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+objectnorder-1, # number of motifs coefficients for each selected motif
            motifs_in_curves_i$starting_coeff_pos-1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
          list_coeff <- lapply(error_str,function(sd_noise){
            coeff[pos_coeff_motifs]=unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
              print(paste(" --- Adding motif", id, "to curve", j,"with noise ",sd_noise))
              object@mot_details[[id]]$weights})) +
              rnorm(length(pos_coeff_motifs),sd=rep(rep(sd_noise,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id],rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1)) # Adding gaussian noise to coefficients
            # For each chosen coefficient add a uniform number and a gaussian noise
            coeff[-pos_coeff_motifs]=object@coeff_min+rbeta(len_i-length(pos_coeff_motifs),0.45,0.45)*(object@coeff_max-object@coeff_min) # All other curve coefficients are randomly generated
            return(coeff)
          })
        }
      }else{
        stop('Wrong \'distrib\'')
      }
      fda_with_error <- NULL
      or_y <- NULL
      if(!is.null(motifs_in_curves_i)) {
        fda_with_error= Map(fda::fd,list_coeff,MoreArgs = list(basisobj=basis_i)) # Fitting curves using such coefficients and basis
        or_y <- Map(mtfd:::generate_curve_vector,fda_with_error)
      }
      return(list(no_error = list(or_coeff = or_coeff,basis = basis_i,no_error_y = or_y_no_error),
                  with_error = list(error_structure = error_str,error_y = or_y)))
    },object@motifs_in_curves,object@norder-1+rep(object@len/object@dist_knots,length.out=object@N),basis,1:object@N,MoreArgs = list(len_motifs),SIMPLIFY=FALSE)
    
  }else{
    fd_curves <- lapply(1:object@N, function(x){
      coeff <- NULL
      len_i <- object@norder-1+object@len/object@dist_knots
      if(object@distribution=='unif'){
        coeff=runif(len_i,min=object@coeff_min,max=object@coeff_max) # If we don't have motifs we sample len_i coefficients uniformely
      }else{
        coeff=object@coeff_min+rbeta(len_i,0.45,0.45)*(object@coeff_max-object@coeff_min) 
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
        t = seq(0, curves[[k]]$no_error$basis$rangeval[2]-1),
        x = curves[[k]]$no_error$no_error_y)
      names(curve_data_no_error) <- c("t","x")
      
      curve_data_error <- NULL
      curve_data_error <- data.frame(
        t = seq(0, curves[[k]]$no_error$basis$rangeval[2]-1),
        x = curves[[k]]$with_error$error_y)
      names(curve_data_error) <- c("t",paste0("x",seq(length(curves[[k]]$with_error$error_y))))
      
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
      p <- lapply(1:length(motif_data),function(j) {
        pic <- ggplot() +
          # Plot the main curve in black
          geom_line(data = curve_data_no_error, aes(x = t, y = x), color = scales::alpha('gray30',0.15), linewidth = 0.5) + 
          # Plot the error curve in black
          geom_line(data = curve_data_error, aes_string(x = "t", y = paste0("x", j)), color = "black", linewidth = 0.5) +
          # Add shaded rectangles for motif positions with transparency
          geom_rect(data = motif_data[[j]], aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = initial_number), alpha = 0.005) +
          # Plot motifs with distinct colors
          geom_line(data = motif_data[[j]], aes(x = t, y = x, color = factor(initial_number), group = motif_id), linewidth = 1.0) + 
          scale_color_manual(values = motif_colors) +
          scale_fill_manual(values = motif_colors) +
          labs(title = paste('Random curve', k), aes_string(x = "t", y = paste0("x", j))) +
          theme_minimal(base_size = 15) +
          guides(color = guide_legend(title = "Motif ID"), fill = guide_legend(title = "Motif ID"))
        return(pic)
      }) 
      Map(print,p)
    }
  }
  dev.off()
})
