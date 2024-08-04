
.transform_to_matrix <- function(input) {
  if (is.vector(input)) {
    if (length(input) == 1) {
      # Se l'input è un singolo numero, trasformalo in una matrice 1x1
      result <- matrix(input, nrow = 1, ncol = 1)
    } else {
      # Se l'input è un vettore, trasformalo in una matrice 1xN
      result <- matrix(input, nrow = 1, ncol = length(input))
    }
  } else if (is.matrix(input)) {
    # Se l'input è già una matrice, non fare nulla
    result <- input
  } else {
    stop("Input must be numeric or a matrix")
  }
  return(result)
}

.transform_list <- function(lst,error_str) {
  # Looping the list
  for (i in seq_along(lst)) {
    # check"with_error" field
    if (!"with_error" %in% names(lst[[i]])) {
      # Create a new sub-list
      lst[[i]] <- list(
        no_error = lst[[i]],
        with_error = list(
          error_structure = mtfd:::.transform_to_matrix(error_str),
          error_y = NULL)
      )
    }
  }
  return(lst)
}

generate_curve_vector <- function(fd_curve, step_by = 1, Lfdobj = 0){
  # from an fd object generate the curve as a vector of len fixed by the fd_obj
  # with step equal to step_by. It generates the derivative Lfdobj
  x <- seq(0, fd_curve$basis$rangeval[2] - 1, by = step_by)
  y <- eval.fd(x, fd_curve, Lfdobj = Lfdobj)  # return the evaluation of the Lfdobj derivative
  return(y)
}

# add an additive noise to the final curve(motif)
add_error_to_motif <- function(or_y, error_str, start_point, end_point,k){
  # --- additive noise
  # noise to be added as a percentage of sd(motif) with mean 0
  if(length(start_point:end_point) - 1 > length(error_str)) {
    
    warning(paste("error structure on row",k,"is long",length(error_str),
                  "but the motif is long",length(start_point:end_point) - 1,". -- Extension"))
    # Last value
    last_value <- tail(error_str, 1)
    
    # Extention 
    error_str <- c(error_str, rep(last_value, length(start_point:end_point) - length(error_str)))
  }else if(length(start_point:end_point) - 1 < length(error_str)) {
    
    warning(paste("error structure on row",k,"is long",length(error_str),
                  "but the motif is long",length(start_point:end_point) - 1,". -- Truncation"))
    #Truncation
    error_str <- error_str[1:length(start_point:end_point)]
  }
  noise_coeff <- rnorm(length(start_point:end_point) - 1,
                       mean = 0,
                       sd = sd(or_y[start_point:end_point])*error_str)
  # add noise
  err_y <- or_y
  err_y[(start_point+1):end_point] <- err_y[(start_point+1):end_point] + noise_coeff
  
  return(err_y)
}

# generate a background curve with length (len) using b-spline with knots at (dist_knots)
# one from the other, order (norder) and weights sampled from (weights). User can 
# decide to add noise (add_noise = TRUE) sampled from a Gaussian distribution
# with mu = 0 and sd = 0.01 (useful for additive noise in motifs)
generate_background_curve <- function(len, dist_knots, norder, weights, add_noise){
  
  # create knots and generate the corresponding b-spline basis for every curve
  breaks <- seq(from = 0, len, by = dist_knots)
  basis  <- create.bspline.basis(norder = norder, breaks = breaks)
  
  if(length(weights) >= basis$nbasis){
    or_coeff <- weights[1:(basis$nbasis)]
  } else {
    print('A longer vector of "weights" is required')
    break
  }
  
  # create fd object
  fd_curve <- fd(or_coeff, basis)
  # create the curve vector with no noise
  or_y <- generate_curve_vector(fd_curve = fd_curve)
  
  if(add_noise == TRUE){
    or_y <- or_y + rnorm(n = length(or_y), mean = 0, sd = 0.01)
  }
  
  # create list with:
  res <- list("or_coeff" = or_coeff,
              "basis" = basis,
              "no_error_y" = or_y # curve as vector with no error
  ) 
  return(res)
}

add_motif <- function(base_curve, mot_pattern, mot_len, dist_knots, mot_order, mot_weights, error_str){
  
  # create knots and generate the corresponding b-spline basis for every curve
  mot_breaks <- seq(from = 0, mot_len, by = dist_knots)
  mot_basis  <- create.bspline.basis(norder = mot_order, breaks = mot_breaks)
  
  if(length(mot_weights) >= mot_basis$nbasis){
    mot_coeff <- mot_weights[1:(mot_basis$nbasis)]
  } else {
    print('A longer vector of "weights" is required')
    break
  }
  
  # add info about motif: id, starting break or point, ending break or point
  base_curve_coeff <- base_curve$or_coeff
  motif_recap <- cbind.data.frame()
  err_y_list <- list()
  add_err_y_list <- list()
  for(i in 1:nrow(mot_pattern)){
    start_break <- mot_pattern[i, "start_break_pos"] # select the first coefficient 
    start_point <- (start_break-1)*dist_knots # select the first knots
    end_break   <- start_break + mot_basis$nbasis - mot_order + 2 # select the last coefficient  DA CAMBIARE
    end_point   <- start_point + mot_len # select the last knots
    motif_recap[i, "motif_id"]   <- mot_pattern[i, "motif_id"]
    motif_recap[i,"start_break_pos"] <- start_break
    motif_recap[i,"start_point"] <- start_point
    motif_recap[i,"end_break_pos"]   <- end_break
    motif_recap[i,"end_point"]   <- end_point
    base_curve_coeff[start_break:end_break] <- mot_coeff
  }
  
  # Errorless curve
  # create fd object
  # Convert base_curve_coeff and basis to fd object
  fd_curve <- fd(base_curve_coeff, base_curve$basis)
  
  # Generate curve vector with no noise
  or_y <- generate_curve_vector(fd_curve = fd_curve)
  
  no_error_res <- list("or_coeff" = base_curve_coeff,
                       "basis" = base_curve$basis,
                       "no_error_y" = or_y # curve as vector with no error
  )
  
  # Error curve
  # add extra error on motif
  err_y_mat <- list()
  error_str <- .transform_to_matrix(error_str)
  for(k in 1:nrow(error_str)) {
    err_y <- or_y
    error_str_k <- error_str[k,]
    for(i in 1:nrow(mot_pattern)){
      start_break <- mot_pattern[i, "start_break_pos"]
      start_point <- (start_break-1)*dist_knots
      end_break   <- start_break + mot_basis$nbasis - mot_order + 2 # mot_basis$nbasis - order + 2 
      end_point   <- start_point + mot_len
      err_y <- mtfd:::add_error_to_motif(err_y, error_str_k, start_point, end_point,k)
    }
    
    # and then smooth it again
    mot_breaks <- seq(from = 0, length(err_y), by = dist_knots)
    basisobj = create.bspline.basis(norder = mot_order, breaks = mot_breaks)
    #plot(basisobj)
    ys = smooth.basis(argvals = 1:(length(err_y)),
                      y = err_y,
                      fdParobj = basisobj) # smooth again given the error on data(pointwise)
    # smoothed again after adding error
    yy <- eval.fd(1:(length(err_y)), ys$fd)
    err_y_mat[[k]] <- yy
    #lines(yy, col='blue')
  }

  error_res <- list("error_structure" = error_str,
                    "error_y" = err_y_mat) 
  res <- list("no_error" = no_error_res,
              "with_error" = error_res)
  # create list with 
  return(res)
}

.generate_coefficients <- function(motif_i, distrib, dist_knots, norder, coeff_min, coeff_max) {
  # Calculate the length of the coefficients vector
  l <- motif_i$len / dist_knots + norder - 1
  if (distrib == "unif") {
    # Generate coefficients from a uniform distribution
    motif_i$weights <- runif(l, min = coeff_min, max = coeff_max)
  } else if (distrib == "beta") {
    # Generate coefficients from a beta distribution and scale to the desired range
    motif_i$weights <- coeff_min + rbeta(l, shape1 = 0.45, shape2 = 0.45) * (coeff_max - coeff_min)
  } else {
    stop("Wrong 'distrib': ", distrib)
  }
  return(motif_i)
}

# Check for overlaps and fitting within each curve
.check_fits <- function(df,min_dist_motifs) {
  df <- df[order(df$start), ]
  # Check for overlaps
  for (i in seq_len(nrow(df) - 1)) {
    if (df$end[i] + min_dist_motifs > df$start[i + 1]) {
      return(FALSE)
    }
  }
  
  # Check if the last motif fits within the curve length
  if (nrow(df) > 0 && df$end[nrow(df)] > 300) {
    return(FALSE)
  }
  
  return(TRUE)
}

.resample <- function(x, ...) x[sample.int(length(x), ...)]



