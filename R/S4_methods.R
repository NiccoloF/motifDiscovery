#' @title Generate Functional Curves with Embedded Motifs
#'
#' @description
#' The `generateCurves` function is designed to create synthetic functional data by embedding predefined motifs into background curves. This is particularly useful for testing and benchmarking motif discovery and clustering algorithms in functional data analysis. By allowing the incorporation of various noise types and controlled motif placements, the function provides a flexible framework for simulating realistic scenarios where motifs may or may not be noisy.
#'
#' The function supports two types of noise addition:
#' \itemize{
#'   \item \strong{Pointwise Noise}: Adds noise directly to the data points of the curves, simulating random fluctuations or measurement errors.
#'   \item \strong{Coefficient Noise}: Perturbs the coefficients of the basis functions used to represent the curves, allowing for smoother variations and controlled distortions.
#' }
#'
#' Additionally, `generateCurves` allows for the specification of vertical shifts to motifs, enabling the simulation of motifs appearing at different baseline levels within the curves. The function ensures that all generated subcurves meet the minimum motif length requirement, maintaining the integrity of the embedded motifs.
#'
#' This function is integral to the `funMoDisco` package's motif simulation capabilities, providing users with the ability to create complex functional datasets tailored to their specific research or testing needs.
#'
#' @param object An S4 object of class `motifSimulation` that has been previously constructed using the `motifSimulationBuilder` function. This object encapsulates all necessary parameters and configurations for curve and motif generation (mandatory).
#' @param noise_type A character string specifying the type of noise to add to the curves. Acceptable values are `'pointwise'` for adding noise directly to data points or `'coeff'` for perturbing the coefficients of the basis functions (mandatory).
#' @param noise_str A list detailing the structure and magnitude of the noise to be added for each motif. 
#'   \itemize{
#'     \item If `noise_type` is `'pointwise'`, `noise_str` should contain vectors or matrices indicating the noise level for each motif.
#'     \item If `noise_type` is `'coeff'`, `noise_str` should include individual values or vectors representing the noise to be applied to the coefficients.
#'   }
#'   This parameter allows fine-grained control over the noise characteristics applied to each motif (mandatory).
#' @param seed_background An integer value setting the seed for the random number generator used in background curve generation. This ensures reproducibility of the background curves. Default is `777`.
#' @param only_der A logical value indicating whether to apply only derivative-based modifications to the motifs (`TRUE`) or to add a vertical shift in addition to derivative modifications (`FALSE`). Setting to `FALSE` introduces vertical shifts to each motif instance, allowing motifs to appear at different baseline levels within the curves. Default is `TRUE`.
#' @param coeff_min_shift A numeric value specifying the minimum vertical shift to be applied to motifs when `only_der` is set to `FALSE`. This parameter controls the lower bound of the vertical displacement of motifs. Default is `-10`.
#' @param coeff_max_shift A numeric value specifying the maximum vertical shift to be applied to motifs when `only_der` is set to `FALSE`. This parameter controls the upper bound of the vertical displacement of motifs. Default is `10`.
#'
#' @return 
#' A list containing the following components:
#' \itemize{
#'   \item \strong{basis}: The basis functions used to represent the curves.
#'   \item \strong{background}: A list containing the coefficients and the background curves without any motifs.
#'   \item \strong{no_noise}: A list containing the coefficients and the background curves with embedded motifs but without added noise.
#'   \item \strong{with_noise}: A list containing the noise structure and the curves with embedded motifs and added noise.
#'   \item \strong{SNR}: A list of Signal-to-Noise Ratio (SNR) metrics calculated for each motif within each curve, useful for assessing the quality of motif embedding.
#' }
#'
#' @examples
#' \donttest{
#' # Example 0: Special case with no motifs
#' mot_len <- 100
#' mot_details <- NULL  # or list()
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details)
#' curves <- funMoDisco::generateCurves(builder)
#'
#' # Example 1: Set the motif position and add pointwise noise
#' # Define motif positions and their respective curves
#' motif_str <- rbind.data.frame(
#'   c(1, 1, 20),
#'   c(2, 1, 2),
#'   c(2, 7, 1),
#'   c(2,17,1)
#' )
#' names(motif_str) <- c("motif_id", "curve", "start_break_pos")
#'
#' # Define motif details
#' mot1 <- list(
#'   "len" = mot_len, 
#'   "coeffs" = NULL, 
#'   "occurrences" = motif_str %>% filter(motif_id == 1)
#' )
#' mot2 <- list(
#'   "len" = mot_len, 
#'   "coeffs" = NULL, 
#'   "occurrences" = motif_str %>% filter(motif_id == 2)
#' )
#' mot_details <- list(mot1, mot2)
#'
#' # Define noise structure for pointwise noise
#' noise_str <- list(
#'   rbind(rep(2, 100), rep(c(rep(0.1, 50), rep(2, 50)), 1)),
#'   rbind(rep(0.0, 100), rep(0.5, 100))
#' )
#'
#' # Build the simulation object
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')
#'
#' # Generate curves with pointwise noise
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str = noise_str)
#'
#' # Example 2: Set the motif position and add coefficient noise
#' # Define noise structure for coefficient noise
#' noise_str <- list(c(0.1, 1.0, 5.0), c(0.0, 0.0, 0.0))
#'
#' # Generate curves with coefficient noise without vertical shifts
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
#'
#' # Example 3: Random motif positions and add pointwise noise
#' mot1 <- list(
#'   "len" = mot_len,
#'   "coeffs" = NULL,
#'   "occurrences" = 5
#' )
#' mot2 <- list(
#'   "len" = mot_len,
#'   "coeffs" = NULL,
#'   "occurrences" = 6
#' )
#' mot_details <- list(mot1, mot2)
#'
#' # Define noise structure for pointwise noise
#' noise_str <- list(
#'   rbind(rep(2, 100)),
#'   rbind(rep(0.5, 100))
#' )
#'
#' # Build the simulation object
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')
#'
#' # Generate curves with pointwise noise and vertical shifts
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str, only_der = FALSE)
#'
#' # Example 4: Random motif positions and add coefficient noise
#' # Define noise structure for coefficient noise
#' noise_str <- list(c(0.1, 5.0, 10.0), c(0.1, 5.0, 10.0))
#'
#' # Generate curves with coefficient noise and vertical shifts
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
#' }
setGeneric("generateCurves", function(object,noise_type = NULL, noise_str = NULL,seed_background = NULL,
                                      only_der = FALSE, coeff_min_shift = -10,coeff_max_shift = 10) {
  standardGeneric("generateCurves")
})


#' @title Generate Functional Curves with Embedded Motifs
#'
#' @description
#' The `generateCurves` function is designed to create synthetic functional data by embedding predefined motifs into background curves. This is particularly useful for testing and benchmarking motif discovery and clustering algorithms in functional data analysis. By allowing the incorporation of various noise types and controlled motif placements, the function provides a flexible framework for simulating realistic scenarios where motifs may or may not be noisy.
#'
#' The function supports two types of noise addition:
#' \itemize{
#'   \item \strong{Pointwise Noise}: Adds noise directly to the data points of the curves, simulating random fluctuations or measurement errors.
#'   \item \strong{Coefficient Noise}: Perturbs the coefficients of the basis functions used to represent the curves, allowing for smoother variations and controlled distortions.
#' }
#'
#' Additionally, `generateCurves` allows for the specification of vertical shifts to motifs, enabling the simulation of motifs appearing at different baseline levels within the curves. The function ensures that all generated subcurves meet the minimum motif length requirement, maintaining the integrity of the embedded motifs.
#'
#' This function is integral to the `funMoDisco` package's motif simulation capabilities, providing users with the ability to create complex functional datasets tailored to their specific research or testing needs.
#'
#' @param object An S4 object of class `motifSimulation` that has been previously constructed using the `motifSimulationBuilder` function. This object encapsulates all necessary parameters and configurations for curve and motif generation (mandatory).
#' @param noise_type A character string specifying the type of noise to add to the curves. Acceptable values are `'pointwise'` for adding noise directly to data points or `'coeff'` for perturbing the coefficients of the basis functions (mandatory).
#' @param noise_str A list detailing the structure and magnitude of the noise to be added for each motif. 
#'   \itemize{
#'     \item If `noise_type` is `'pointwise'`, `noise_str` should contain vectors or matrices indicating the noise level for each motif.
#'     \item If `noise_type` is `'coeff'`, `noise_str` should include individual values or vectors representing the noise to be applied to the coefficients.
#'   }
#'   This parameter allows fine-grained control over the noise characteristics applied to each motif (mandatory).
#' @param seed_background An integer value setting the seed for the random number generator used in background curve generation. This ensures reproducibility of the background curves. Default is `777`.
#' @param only_der A logical value indicating whether to apply only derivative-based modifications to the motifs (`TRUE`) or to add a vertical shift in addition to derivative modifications (`FALSE`). Setting to `TRUE` introduces vertical shifts to each motif instance, allowing motifs to appear at different baseline levels within the curves. Default is `FALSE`.
#' @param coeff_min_shift A numeric value specifying the minimum vertical shift to be applied to motifs when `only_der` is set to `TRUE`. This parameter controls the lower bound of the vertical displacement of motifs. Default is `-10`.
#' @param coeff_max_shift A numeric value specifying the maximum vertical shift to be applied to motifs when `only_der` is set to `TRUE`. This parameter controls the upper bound of the vertical displacement of motifs. Default is `10`.
#'
#' @return 
#' A list containing the following components:
#' \itemize{
#'   \item \strong{basis}: The basis functions used to represent the curves.
#'   \item \strong{background}: A list containing the coefficients and the background curves without any motifs.
#'   \item \strong{no_noise}: A list containing the coefficients and the background curves with embedded motifs but without added noise.
#'   \item \strong{with_noise}: A list containing the noise structure and the curves with embedded motifs and added noise.
#'   \item \strong{SNR}: A list of Signal-to-Noise Ratio (SNR) metrics calculated for each motif within each curve, useful for assessing the quality of motif embedding.
#' }
#'
#' @export 
#' 
#' @importFrom purrr is_empty
#' @examples
#' \donttest{
#' # Example 0: Special case with no motifs
#' mot_len <- 100
#' mot_details <- NULL  # or list()
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details)
#' curves <- funMoDisco::generateCurves(builder)
#'
#' # Example 1: Set the motif position and add pointwise noise
#' # Define motif positions and their respective curves
#' motif_str <- rbind.data.frame(
#'   c(1, 1, 20),
#'   c(2, 1, 2),
#'   c(2, 7, 1),
#'   c(2,17,1)
#' )
#' names(motif_str) <- c("motif_id", "curve", "start_break_pos")
#'
#' # Define motif details
#' mot1 <- list(
#'   "len" = mot_len, 
#'   "coeffs" = NULL, 
#'   "occurrences" = motif_str %>% filter(motif_id == 1)
#' )
#' mot2 <- list(
#'   "len" = mot_len, 
#'   "coeffs" = NULL, 
#'   "occurrences" = motif_str %>% filter(motif_id == 2)
#' )
#' mot_details <- list(mot1, mot2)
#'
#' # Define noise structure for pointwise noise
#' noise_str <- list(
#'   rbind(rep(2, 100), rep(c(rep(0.1, 50), rep(2, 50)), 1)),
#'   rbind(rep(0.0, 100), rep(0.5, 100))
#' )
#'
#' # Build the simulation object
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')
#'
#' # Generate curves with pointwise noise
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str = noise_str)
#'
#' # Example 2: Set the motif position and add coefficient noise
#' # Define noise structure for coefficient noise
#' noise_str <- list(c(0.1, 1.0, 5.0), c(0.0, 0.0, 0.0))
#'
#' # Generate curves with coefficient noise without vertical shifts
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
#'
#' # Example 3: Random motif positions and add pointwise noise
#' mot1 <- list(
#'   "len" = mot_len,
#'   "coeffs" = NULL,
#'   "occurrences" = 5
#' )
#' mot2 <- list(
#'   "len" = mot_len,
#'   "coeffs" = NULL,
#'   "occurrences" = 6
#' )
#' mot_details <- list(mot1, mot2)
#'
#' # Define noise structure for pointwise noise
#' noise_str <- list(
#'   rbind(rep(2, 100)),
#'   rbind(rep(0.5, 100))
#' )
#'
#' # Build the simulation object
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')
#'
#' # Generate curves with pointwise noise and vertical shifts
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str, only_der = FALSE)
#'
#' # Example 4: Random motif positions and add coefficient noise
#' # Define noise structure for coefficient noise
#' noise_str <- list(c(0.1, 5.0, 10.0), c(0.1, 5.0, 10.0))
#'
#' # Generate curves with coefficient noise and vertical shifts
#' curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
#' }
setMethod("generateCurves", "motifSimulation", function(object,
                                                        noise_type = NULL,
                                                        noise_str = NULL,
                                                        seed_background = NULL,
                                                        only_der = FALSE,
                                                        coeff_min_shift = -10,
                                                        coeff_max_shift = 10){
  
  breaks = lapply(rep(object@len, length.out = object@N), seq, from = 0, by = object@dist_knots) # generate N lists with equally spaced nodes 'dist_knots' from 0 to 'len'
  basis = lapply(breaks, function(breaks_i) create.bspline.basis(norder = object@norder, breaks = breaks_i)) # nbasis = norder + nbreaks - 2 = norder + interior_knots,
  len_motifs <- unlist(lapply(object@mot_details, function(x){x$len}))
  
  # if the seed is user defined (not null) then set the seed
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
    stop("\'noise_type\' must be choosen between \'coeff\' and \'pointwise\'")
  } 
  if(is.null(noise_type)) {
    noise_type <- "coeff" # to generate background curves without noise
  }
  
  if(noise_type == 'coeff') {
    
    tryCatch({
      # Attempt to coerce the input to a numeric vector
      noise_str <- lapply(noise_str, function(error){return(as.numeric(error))})
    }, error = function(e) {
      stop("\'noise_str'\ cannot be vectorized.")
    })
    
    if(length(noise_str) != length(object@mot_details)) {
      stop("\'noise_str\' must have the same length of \'mot_details\'") 
    }
    
    fd_curves <- mapply(function(motifs_in_curves_i,len_i,basis_i,j,len_motifs){
      list_coeff_error <- NULL
      if(is.numeric(object@distribution)) {
        or_coeff = sample(object@distribution,size = len_i,replace = TRUE)
      } else if(object@distribution == 'unif'){
        or_coeff = runif(len_i, min = object@coeff_min, max = object@coeff_max) # coefficients of the curve = degree of freedom are initialized uniformly
      } else if(object@distribution == 'beta') { 
        or_coeff = object@coeff_min + rbeta(len_i,0.45,0.45)*(object@coeff_max - object@coeff_min) # coefficients of the curve = degree of freedom are initialized uniformly
      } else {
        stop('Wrong \'distrib\'')
      }
      coeff <- or_coeff
      fda_no_error = fda::fd(coef = coeff, basisobj = basis_i)
      or_y_no_error <- generate_curve_vector(fda_no_error)
      fda_with_error <- NULL
      or_y <- NULL
      shifted_coeff <- NULL
      #set.seed(seed_motif)
      if(!is.null(motifs_in_curves_i)) {
        pos_coeff_motifs = unlist(mapply(
          function(a,b) seq(a) + b,
          rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
          motifs_in_curves_i$starting_coeff_pos - 1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
        shifted_coeff <- rep(runif(length(motifs_in_curves_i$motif_id), min = coeff_min_shift, max = coeff_max_shift),rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1) 
        
        list_coeff_error <- lapply(1:length(noise_str[[1]]),function(z) {
          if(only_der == FALSE) {
            coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
              sd_noise <- noise_str[[id]][z]
              message(paste("\033[34m --- Adding motif", id, "to curve", j, "with noise", sd_noise, "\033[0m"))
              object@mot_details[[id]]$coeffs + rnorm(length(object@mot_details[[id]]$coeffs),sd=sd_noise)}))
          } else { 
            coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
              sd_noise <- noise_str[[id]][z]
              message(paste("\033[34m --- Adding motif", id, "to curve", j, "with noise", sd_noise, "\033[0m"))
              object@mot_details[[id]]$coeffs + rnorm(length(object@mot_details[[id]]$coeffs),sd=sd_noise)})) + shifted_coeff
          }
          # For each chosen coefficient add a uniform number and a gaussian noise
          return(coeff)
        })
        if(only_der == FALSE) {
          coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
            object@mot_details[[id]]$coeffs}))
        } else {
          coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
            object@mot_details[[id]]$coeffs})) + shifted_coeff
        }
        fda_with_motif <- fda::fd(coef = coeff, basisobj = basis_i)
        or_y_motif <- generate_curve_vector(fda_with_motif)
        
        fda_with_error <- Map(fda::fd,list_coeff_error,MoreArgs = list(basisobj=basis_i)) # Fitting curves using such coefficients and basis
        or_y <- Map(generate_curve_vector,fda_with_error)
        num_motifs <- length(motifs_in_curves_i$motif_id)
        SNR <- vector("list",length(noise_str[[1]]))
        for (k in seq_along(noise_str[[1]])) {
          SNR_num <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
          SNR_den <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
          
          # Loop over each motif
          for (n in seq_along(motifs_in_curves_i$motif_id)) {
            start_break <- motifs_in_curves_i$starting_coeff_pos[n]
            start_point <- (start_break - 1) * object@dist_knots
            end_point   <- start_point + len_motifs[motifs_in_curves_i$motif_id[n]]
            
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
    
  }else if(noise_type == 'pointwise'){
    
    if(length(noise_str) != length(object@mot_details)) {
      stop("\'noise_str\' must have the same length of \'mot_details\'") 
    }
    
    fd_curves <- lapply(1:object@N, function(x){
      coeff <- NULL
      len_i <- object@norder - 1 + object@len/object@dist_knots
      if(is.numeric(object@distribution)) {
        coeff = sample(object@distribution, size = len_i, replace = TRUE)
      }
      else if(object@distribution=='unif'){
        coeff = runif(len_i, min = object@coeff_min, max = object@coeff_max) # If we don't have motifs we sample len_i coefficients uniformely
      }else if(object@distribution=='beta'){
        coeff = object@coeff_min + rbeta(len_i, 0.45, 0.45)*(object@coeff_max - object@coeff_min) 
      }else {
        stop('Wrong \'distrib\'')
      }
      generate_background_curve(object@len, object@dist_knots, object@norder, coeff, add_noise = TRUE)
    })
    
    #set.seed(seed_motif)
    curve_ids <- unique(unlist(lapply(object@mot_details,function(mot){mot$occurrences$curve})))
    for(j in curve_ids){
      message(paste("\033[34m --- Adding motifs to curve", j, "\033[0m"))
      temp_curve <- fd_curves[[j]]
      temp_pattern <- do.call(rbind,lapply(object@mot_details,function(mot_details){(mot_details$occurrences %>% filter(curve == j)) %>% dplyr::select(motif_id, start_break_pos)}))
      temp_len <- do.call(rbind,lapply(unique(apply(temp_pattern,1,function(row){row[1]})),function(k){data.frame("motif_id" = k,"len" = object@mot_details[[k]]$len)}))
      temp_weights <- lapply(as.character(unique(apply(temp_pattern, 1, function(row) { row[1] }))), 
                             function(k) { 
                               object@mot_details[[as.numeric(k)]]$coeffs
                             })
      names(temp_weights) <- as.character(unique(apply(temp_pattern, 1, function(row) { row[1] })))
      fd_curves[[j]] <- add_motif(base_curve  = temp_curve,
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
    fd_curves <- .transform_list(fd_curves,noise_str)
  } else {
    stop("\'noise_type\' must be choosen between \'coeff\' and \'pointwise\'")
  }
  return(fd_curves = fd_curves)
})

# PREVIOUS VERSION
# setMethod("generateCurves", "motifSimulation", function(object, noise_type = NULL, noise_str = NULL, seed_background = NULL, seed_motif = NULL,
#                                                         only_der = FALSE, coeff_min_shift = -10, coeff_max_shift = 10) {
#   breaks = lapply(rep(object@len, length.out = object@N), seq, from = 0, by = object@dist_knots) # generate N lists with equally spaced nodes 'dist_knots' from 0 to 'len'
#   basis = lapply(breaks, function(breaks_i) create.bspline.basis(norder = object@norder, breaks = breaks_i)) # nbasis = norder + nbreaks - 2 = norder + interior_knots,
#   len_motifs <- unlist(lapply(object@mot_details, function(x){x$len}))
#   
#   # if the seed is user defined (not null) then set the seed
#   if(!is.null(seed_background)){
#     set.seed(seed_background)
#   }
# 
#   # loop for each curve
#   fd_curves <- NULL
#   if(is.null(noise_type) && !is.null(noise_str)) {
#     stop("\'noise_str'\ must be specified because \'noise_type\' is not null")
#   }
#   if(is.null(noise_str) && !is.null(noise_type)) {
#     stop("\'noise_type'\ must be specified because \'noise_str\' is not null")
#   }
#   if((is.null(noise_type) && !is_empty(object@motifs_in_curves))) {
#     stop("\'noise_type\' must be choosen between \'coeff\' and \'pointiwse\'")
#   } 
#   if(is.null(noise_type)) {
#     noise_type <- "pointwise" # to generate background curves without noise
#   }
#   if(noise_type == 'coeff') {
#     tryCatch({
#       # Attempt to coerce the input to a numeric vector
#       noise_str <- lapply(noise_str, function(error){return(as.numeric(error))})
#     }, error = function(e) {
#       stop("\'noise_str'\ cannot be vectorized.")
#     })
#     if(length(noise_str) != length(object@mot_details)) {
#      stop("\'noise_str\' must have the same length of \'mot_details\'") 
#     }
#     fd_curves <- mapply(function(motifs_in_curves_i,len_i,basis_i,j,len_motifs){
#       list_coeff_error <- NULL
#       if(is.numeric(object@distribution)) {
#           or_coeff = sample(object@distribution,size = len_i,replace = TRUE)
#       } else if(object@distribution == 'unif'){
#           or_coeff = runif(len_i, min = object@coeff_min, max = object@coeff_max) # coefficients of the curve = degree of freedom are initialized uniformly
#       } else if(object@distribution == 'beta') { 
#           or_coeff = object@coeff_min + rbeta(len_i,0.45,0.45)*(object@coeff_max - object@coeff_min) # coefficients of the curve = degree of freedom are initialized uniformly
#       } else {
#         stop('Wrong \'distrib\'')
#       }
#       coeff <- or_coeff
#       fda_no_error = fda::fd(coef = coeff, basisobj = basis_i)
#       or_y_no_error <- generate_curve_vector(fda_no_error)
#       fda_with_error <- NULL
#       or_y <- NULL
#       shifted_coeff <- NULL
#       #set.seed(seed_motif)
#       if(!is.null(motifs_in_curves_i)) {
#         pos_coeff_motifs = unlist(mapply(
#           function(a,b) seq(a) + b,
#           rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1, # number of motifs coefficients for each selected motif
#           motifs_in_curves_i$starting_coeff_pos - 1,SIMPLIFY=FALSE)) # Calculating the position of the coefficients of the motifs within the coefficients of the curves
#         shifted_coeff <- rep(runif(length(motifs_in_curves_i$motif_id), min = coeff_min_shift, max = coeff_max_shift),rep(len_motifs,length.out=length(object@mot_details))[motifs_in_curves_i$motif_id]/object@dist_knots+object@norder-1) 
#       
#         list_coeff_error <- lapply(1:length(noise_str[[1]]),function(z) {
#           if(only_der == FALSE) {
#           coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
#             sd_noise <- noise_str[[id]][z]
#             message(paste("\033[34m --- Adding motif", id, "to curve", j, "with noise", sd_noise, "\033[0m"))
#             object@mot_details[[id]]$coeffs + rnorm(length(object@mot_details[[id]]$coeffs),sd=sd_noise)}))
#           } else { 
#             coeff[pos_coeff_motifs] = unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
#               sd_noise <- noise_str[[id]][z]
#               message(paste("\033[34m --- Adding motif", id, "to curve", j, "with noise", sd_noise, "\033[0m"))
#               object@mot_details[[id]]$coeffs + rnorm(length(object@mot_details[[id]]$coeffs),sd=sd_noise)})) + shifted_coeff
#           }
#           # For each chosen coefficient add a uniform number and a gaussian noise
#           return(coeff)
#         })
#       if(only_der == FALSE) {
#         coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
#                                   object@mot_details[[id]]$coeffs}))
#       } else {
#         coeff[pos_coeff_motifs] <- unlist(lapply(motifs_in_curves_i$motif_id, function(id) {
#                                     object@mot_details[[id]]$coeffs})) + shifted_coeff
#       }
#       fda_with_motif <- fda::fd(coef = coeff, basisobj = basis_i)
#       or_y_motif <- generate_curve_vector(fda_with_motif)
#       
#       fda_with_error <- Map(fda::fd,list_coeff_error,MoreArgs = list(basisobj=basis_i)) # Fitting curves using such coefficients and basis
#       or_y <- Map(generate_curve_vector,fda_with_error)
#       num_motifs <- length(motifs_in_curves_i$motif_id)
#       SNR <- vector("list",length(noise_str[[1]]))
#       for (k in seq_along(noise_str[[1]])) {
#         SNR_num <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
#         SNR_den <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
#       
#         # Loop over each motif
#         for (n in seq_along(motifs_in_curves_i$motif_id)) {
#           start_break <- motifs_in_curves_i$starting_coeff_pos[n]
#           start_point <- (start_break - 1) * object@dist_knots
#           end_point   <- start_point + len_motifs[motifs_in_curves_i$motif_id[n]]
# 
#           # Calculate variance for SNR numerator and denominator
#           SNR_num <- rbind(SNR_num,data.frame(xmin =start_point,xmax = end_point,SNR = var(or_y_motif[start_point:end_point])))
#           SNR_den <- rbind(SNR_den,data.frame(xmin =start_point,xmax = end_point,SNR = var(or_y[[k]][start_point:end_point] - or_y_motif[start_point:end_point])))
#         }
#         # Calculate SNR in decibels for the k-th error
#         SNR[[k]] <- SNR_num
#         SNR[[k]]$SNR <- 10 * log10(SNR_num$SNR / SNR_den$SNR) #transform SNR in decibel
#       }
#       return(list(basis = basis_i,
#                   background = list(or_coeff = or_coeff,no_error_y = or_y_no_error),
#                   no_noise = list(or_coeff = coeff,motif_y = or_y_motif),
#                   with_noise = list(noise_structure = noise_str,noise_y = or_y),
#                   SNR = SNR))
#       }
#         
#       return(list(basis = basis_i,
#                   background = list(or_coeff = or_coeff,no_error_y = or_y_no_error)))
#     },object@motifs_in_curves,object@norder-1+rep(object@len/object@dist_knots,length.out=object@N),basis,1:object@N,MoreArgs = list(len_motifs),SIMPLIFY=FALSE)
#     
#   }else if(noise_type == 'pointwise'){
#     if(length(noise_str) != length(object@mot_details)) {
#       stop("\'noise_str\' must have the same length of \'mot_details\'") 
#     }
#     fd_curves <- lapply(1:object@N, function(x){
#       coeff <- NULL
#       len_i <- object@norder - 1 + object@len/object@dist_knots
#       if(is.numeric(object@distribution)) {
#         coeff = sample(object@distribution, size = len_i, replace = TRUE)
#       }
#       else if(object@distribution=='unif'){
#         coeff = runif(len_i, min = object@coeff_min, max = object@coeff_max) # If we don't have motifs we sample len_i coefficients uniformely
#       }else if(object@distribution=='beta'){
#         coeff = object@coeff_min + rbeta(len_i, 0.45, 0.45)*(object@coeff_max - object@coeff_min) 
#       }else {
#         stop('Wrong \'distrib\'')
#       }
#       generate_background_curve(object@len, object@dist_knots, object@norder, coeff, add_noise = TRUE)
#     })
#     
#     #set.seed(seed_motif)
#     curve_ids <- unique(unlist(lapply(object@mot_details,function(mot){mot$occurrences$curve})))
#     for(j in curve_ids){
#       message(paste("\033[34m --- Adding motifs to curve", j, "\033[0m"))
#       temp_curve <- fd_curves[[j]]
#       temp_pattern <- do.call(rbind,lapply(object@mot_details,function(mot_details){(mot_details$occurrences %>% filter(curve == j)) %>% dplyr::select(motif_id, start_break_pos)}))
#       temp_len <- do.call(rbind,lapply(unique(apply(temp_pattern,1,function(row){row[1]})),function(k){data.frame("motif_id" = k,"len" = object@mot_details[[k]]$len)}))
#       temp_weights <- lapply(as.character(unique(apply(temp_pattern, 1, function(row) { row[1] }))), 
#                              function(k) { 
#                                object@mot_details[[as.numeric(k)]]$coeffs
#                              })
#       names(temp_weights) <- as.character(unique(apply(temp_pattern, 1, function(row) { row[1] })))
#       fd_curves[[j]] <- add_motif(base_curve  = temp_curve,
#                                   mot_pattern = temp_pattern,
#                                   mot_len     = temp_len,
#                                   dist_knots  = object@dist_knots,
#                                   mot_order   = object@norder,
#                                   mot_weights = temp_weights,
#                                   noise_str   = noise_str,
#                                   only_der = only_der,
#                                   coeff_min_shift = coeff_min_shift,
#                                   coeff_max_shift = coeff_max_shift)
#     }
#     fd_curves <- .transform_list(fd_curves,noise_str)
#   } else {
#     stop("\'noise_type\' must be choosen between \'coeff\' and \'pointwise\'")
#   }
#   return(fd_curves = fd_curves)
# })

#' @title Plot Embedded Motifs in Functional Curves
#'
#' @description
#' The `plot_motifs` function visualizes the results generated by the `generateCurves` function. It provides comprehensive plots of functional curves with embedded motifs, allowing users to inspect the placement and characteristics of each motif within the curves. This visualization is crucial for validating the correctness of motif embedding and for gaining insights into the distribution and variability of motifs across different curves.
#'
#' The function supports saving the generated plots to a specified directory, facilitating the creation of reports or the sharing of visual results. By plotting motifs in distinct colors or styles, `plot_motifs` ensures that overlapping motifs and their respective curves remain distinguishable, enhancing the clarity and interpretability of the visualizations.
#'
#' This plotting utility is an essential tool within the `funMoDisco` package, aiding users in the exploratory analysis of simulated functional data and the assessment of motif detection algorithms.
#'
#' @param object An S4 object of class `motifSimulation` that has been previously constructed using the `motifSimulationBuilder` function. This object contains all necessary parameters and configurations used during curve and motif generation (mandatory).
#' @param name Name of the output file.
#' @param curves The output list from the `generateCurves` function, containing the generated functional curves with embedded motifs. This parameter provides the data to be visualized (mandatory).
#' @param path A character string specifying the directory path where the generated plots will be saved. The function will save the plots in this directory, allowing for organized storage and easy access to visual results (mandatory).
#'
#' @return 
#' The function does not return any value but generates and saves plots of the functional curves with embedded motifs in the specified directory. Each plot visually represents the motifs within the curves, aiding in the qualitative assessment of motif embedding.
#'
#' @export
#' 
#' @importFrom stats setNames
#' @importFrom ggtext element_markdown
#'
#' @examples
#' \donttest{
#' # Example: Plotting motifs in generated curves
#' # Assume 'builder' has been created and 'curves' have been generated using generateCurves
#' 
#' mot_details <- NULL
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details = mot_details)
#' curves <- funMoDisco::generateCurves(builder)
#'
#' # Specify the directory to save plots
#' plots_name <- "plots_1"  
#'
#' # Generate and save the plots
#' funMoDisco::plot_motifs(builder, curves, plots_name, path = tempdir())
#' }
setGeneric("plot_motifs", function(object,curves,name,path) 
  standardGeneric("plot_motifs")
)

#' @title Plot Embedded Motifs in Functional Curves
#'
#' @description
#' The `plot_motifs` function visualizes the results generated by the `generateCurves` function. It provides comprehensive plots of functional curves with embedded motifs, allowing users to inspect the placement and characteristics of each motif within the curves. This visualization is crucial for validating the correctness of motif embedding and for gaining insights into the distribution and variability of motifs across different curves.
#'
#' The function supports saving the generated plots to a specified directory, facilitating the creation of reports or the sharing of visual results. By plotting motifs in distinct colors or styles, `plot_motifs` ensures that overlapping motifs and their respective curves remain distinguishable, enhancing the clarity and interpretability of the visualizations.
#'
#' This plotting utility is an essential tool within the `funMoDisco` package, aiding users in the exploratory analysis of simulated functional data and the assessment of motif detection algorithms.
#'
#' @param object An S4 object of class `motifSimulation` that has been previously constructed using the `motifSimulationBuilder` function. This object contains all necessary parameters and configurations used during curve and motif generation (mandatory).
#' @param curves The output list from the `generateCurves` function, containing the generated functional curves with embedded motifs. This parameter provides the data to be visualized (mandatory).
#' @param name Name of the output file.
#' @param path A character string specifying the directory path where the generated plots will be saved. The function will save the plots in this directory, allowing for organized storage and easy access to visual results (mandatory).
#'
#' @return 
#' The function does not return any value but generates and saves plots of the functional curves with embedded motifs in the specified directory. Each plot visually represents the motifs within the curves, aiding in the qualitative assessment of motif embedding.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Example: Plotting motifs in generated curves
#' # Assume 'builder' has been created and 'curves' have been generated using generateCurves
#' mot_details <- NULL
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details)
#' curves <- funMoDisco::generateCurves(builder)
#'
#' # Specify the directory to save plots
#' plots_name <- "plots_1"  
#'
#' # Generate and save the plots
#' funMoDisco::plot_motifs(builder, curves, plots_name, path = tempdir())
#' }
setMethod("plot_motifs","motifSimulation",
          function(object, curves, name, path) {
            
  # Ensure the path is absolute
  if (!grepl("^(/|[A-Za-z]:\\\\)", path)) {  
    path <- file.path(getwd(), path)  # Convert to absolute path
  }
  
  # Create the directory if it does not exist
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # Ensure the filename has a .pdf extension
  if (!grepl("\\.pdf$", name)) {     
    name <- paste0(name, ".pdf")   
  }
  
  # Construct the output file path
  output_file <- file.path(path, name)
            
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
          title = paste0('<b><span style="color:#0073C2;"> Curve ', k, '</span></b>'),
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
                                                          initial_number = stringr::str_extract(as.character(id_motif), "^[^_]+"),
                                                          xmin = (pos_motif - 1) * object@dist_knots,
                                                          xmax = (pos_motif - 1) * object@dist_knots + object@mot_details[[id_motif]]$len)}))
      }, object@motifs_in_curves[[k]]$motif_id, object@motifs_in_curves[[k]]$starting_coeff_pos, seq_along(object@motifs_in_curves[[k]]$motif_id), SIMPLIFY = FALSE)
      
      motif_colors <- c( "1" = "red", "2" = "blue", "3" = "darkgreen", "4" = "orange",
                         "5" = "purple", "6" = "cyan", "7" = "magenta", "8" = "brown",
                         "9" = "pink", "10" = "grey")
      motif_colors <- rep(motif_colors,length.out = length(object@mot_details))
      if( length(object@mot_details) > 10 )
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
          geom_line(data = curve_data_no_error, aes(x = t, y = z, color = 'zero_noise_curve'), linewidth = 0.5) +
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
                       'zero_noise_curve' = 'gold', 
                       motif_colors),
            labels = c('background_curve' = 'background_curve', 
                       'zero_noise_curve' = 'zero_noise_curve', 
                       setNames(motif_labels, unique(motif_data[[j]]$initial_number)))
          ) +
          scale_fill_manual(
            values = motif_colors, 
            labels = setNames(motif_labels, unique(motif_data[[j]]$initial_number))
          ) +
          # Title
          labs(
            title = paste0('<b><span style="color:#0073C2;">Curve ', k, ' - type_error ', 
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
    curves_data_noise <- vector("list",max_dataframes)
    for(error_n in 1:max_dataframes) {
      temp <- list()
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
                
                curve_y_noise <- curves[[curve_k]]$with_noise$noise_y[[error_n]][start:(start + object@mot_details[[motif_id]]$len - 1)]
                # Creazione di una sequenza x per il plotting
                x <- seq_along(curve_y_noise)
                
                # Aggiungi i dati alla lista
                temp[[length(temp) + 1]] <- data.frame(x = x, y = curve_y_noise, id = motif_id, instance = instance_id) 
              }
            }
          }
        } 
      }
      curves_data_noise[[error_n]] <- temp
    }
  
    # Combina i dati in un unico dataframe
    curves_df <- lapply(curves_data_noise,bind_rows)
    # Ottieni la lista unica degli ID
    unique_ids <- unique(curves_df[[1]]$id)
    
    # Compute the mean motif for each error
    motif_y_means <- vector("list", length(unique_ids))
    index <- 0
    # Iterate over each unique ID
    for (id in unique_ids) {
      index <- index + 1
      # Iterate through each sublist in builder@mot_details for the current id
      for (i in object@mot_details[[index]]$occurrences$curve %>% unique()) {
        curve <- curves[[i]]
        motif_in_curve_i <- object@motifs_in_curves[[i]]
        for(z in 1:length(motif_in_curve_i$motif_id))
        {
          if(motif_in_curve_i$motif_id[z] == id) {
            start <- (motif_in_curve_i$starting_coeff_pos[z] - 1) * object@dist_knots + 1
            # Compute pointwise means across all curves in the subcurve list
            
            motif_y_means[[index]] <- curve$no_noise$motif_y[start:(start + object@mot_details[[index]]$len - 1)]
            
            # MODIFIED BY ME: this code was computing the mean of the noisy occurrences
            # if (is.null(motif_y_means[[index]])) {
            #   motif_y_means[[index]] <- curve$no_noise$motif_y[start:(start + object@mot_details[[index]]$len - 1)]
            # } else {
            #   motif_y_means[[index]] <- motif_y_means[[index]] + curve$no_noise$motif_y[start:(start + object@mot_details[[index]]$len - 1)]
            # }
          }
        }
      }
      #motif_y_means[[index]] <- data.frame(x = seq_along(motif_y_means[[index]]), y = motif_y_means[[index]] / length(object@mot_details[[index]]$occurrences$curve %>% unique()))
      motif_y_means[[index]] <- data.frame(x = seq_along(motif_y_means[[index]]), y = motif_y_means[[index]] %>% unique())
      names(motif_y_means)[index] <- id
    }
    
    # Loop per plottare ogni ID separatamente su pagine diverse
    for (id in unique_ids) {
      for(error_n in 1:max_dataframes) {
        # Filtra i dati per il singolo ID
        plot_data <- curves_df[[error_n]] %>% filter(id == !!id)
        
        p <- ggplot() + 
          geom_line(data = plot_data, aes(x = x, y = y, color = instance), size = 1, linetype = "longdash") + 
          geom_line(data = motif_y_means[[as.character(id)]], aes(x = x, y = y, color = "black"), size = 2, linetype = "solid") + 
          scale_color_manual(
            values = c("motif_y_means" = "black", setNames(rainbow(length(unique(plot_data$instance))), unique(plot_data$instance))),
            labels = c("motif_y_means" = paste0("motif ", id), unique(plot_data$instance))
          ) + 
          labs(
            title = paste0('<b><span style="color:#0073C2;">Motif ', id,'</span></b>'),
            subtitle = paste0("Mean of the no noise motif <b>", id, "</b> with type error <b>", error_n, "</b>"),
            x = "x",
            y = "t",
            color = "Mot. occurrences"
          ) + 
          theme_minimal(base_size = 15) + 
          theme(
            plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 10)),
            plot.subtitle = element_markdown(size = 16, hjust = 0.5, margin = margin(b = 10)),  # Styled subtitle
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
            color = guide_legend(ncol = 1, byrow = TRUE, title = "Mot. occurrences"),
            fill = "none"
          )
        print(p)
      }
    }
  }
  dev.off()
})

#' @title to_motifDiscovery
#' @description Transforms the result of generateCurves into a format suitable for discoverMotifs
#' @param curves A list coming from the generateCurves function.
#' @return A list containing all curves formatted to be suitable for input into the discoverMotifs function.
#' @export
#' @examples
#' \donttest{
#' N <- 20 # number of curves
#' len <- 300 # length of the curves
#' mot_details <- NULL # or list()
#'
#' builder <- motifSimulationBuilder(N = N, len = len, mot_details = NULL)
#' curves <- funMoDisco::generateCurves(builder)
#' formatted_curves <- to_motifDiscovery(curves)
#' }
setGeneric("to_motifDiscovery", function(curves) 
  standardGeneric("to_motifDiscovery")
)

#' @title to_motifDiscovery
#' @description Transforms the result of generateCurves into a format suitable for discoverMotifs
#' @param curves A list coming from the generateCurves function.
#' @return A list containing all curves formatted to be suitable for input into the discoverMotifs function.
#' @export
#' @examples
#' \donttest{
#' builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details = NULL)
#' curves <- funMoDisco::generateCurves(builder)
#' formatted_curves <- to_motifDiscovery(curves)
#' }
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
