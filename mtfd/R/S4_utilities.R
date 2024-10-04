
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

.transform_list <- function(lst,noise_str) {
  # Looping the list
  for (i in seq_along(lst)) {
    # check"with_noise" field
    if (!"with_noise" %in% names(lst[[i]])) {
      # Create a new sub-list
      background <- list(or_coeff = lst[[i]]$or_coeff,    
                         no_error_y = lst[[i]]$no_error_y)
      lst[[i]] <- list(
        basis = lst[[i]]$basis,
        background = background)
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
add_error_to_motif <- function(or_y, noise_str, start_point, end_point,k){
  err_y <- or_y
  # --- additive noise
  # noise to be added as a percentage of sd(motif) with mean 0
  error_str_temp <- noise_str
  for(n in seq_along(start_point)) {
    if(length(start_point[n]:end_point[n]) - 1 > length(noise_str)) {
      
      warning(paste("error structure on row",k,"is long",length(noise_str),
                    "but the motif is long",length(start_point[n]:end_point[n]) - 1,". -- Extension"))
      # Last value
      last_value <- tail(noise_str, 1)
      
      # Extention 
      error_str_temp <- c(noise_str, rep(last_value, length(start_point[n]:end_point[n]) - length(noise_str)))
    }else if(length(start_point[n]:end_point[n]) - 1 < length(noise_str)) {
      
      warning(paste("error structure on row",k,"is long",length(noise_str),
                    "but the motif is long",length(start_point[n]:end_point[n]) - 1,". -- Truncation"))
      #Truncation
      error_str_temp <- noise_str[1:length(start_point[n]:end_point[n])]
    }
    
    noise_coeff <- rnorm(length(start_point[n]:end_point[n]) - 1,
                         mean = 0,
                         sd = sd(or_y[start_point[n]:end_point[n]])*error_str_temp)
  # add noise
    err_y[(start_point[n]+1):end_point[n]] <- err_y[(start_point[n]+1):end_point[n]] + noise_coeff
  }
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
  res <- list("basis" = basis,
              "or_coeff" = or_coeff,
              "no_error_y" = or_y # curve as vector with no error
  ) 
  return(res)
}

add_motif <- function(base_curve, mot_pattern, mot_len, dist_knots, mot_order, mot_weights, noise_str,
                      only_der,coeff_min_shift,coeff_max_shift){
  # create knots and generate the corresponding b-spline basis for every curve
  max_len <- max(mot_len[2])
  mot_breaks <- seq(from = 0,max_len,by = dist_knots)
  mot_basis  <- create.bspline.basis(norder = mot_order, breaks = mot_breaks)
  
  mot_coeff <- vector("list")
  for(j in names(mot_weights)) {
    if(length(mot_weights[[j]]) >= floor(mot_len[mot_len[,1]==as.numeric(j),2] / dist_knots) + 1 + mot_order - 2) { #nbreaks + norder - 2
      mot_coeff[[j]] <- mot_weights[[j]][1:(floor(mot_len[mot_len[,1]==as.numeric(j),2] / dist_knots) + 1 + mot_order - 2)]
    }
    else('A longer vector of "weights" is required')
  }
  # add info about motif: id, starting break or point, ending break or point
  base_curve_coeff <- base_curve$or_coeff
  base_y <- generate_curve_vector(fd_curve = fd(base_curve_coeff, base_curve$basis))
  for(i in names(mot_weights)){
    start_break <- mot_pattern[mot_pattern[,1]==as.numeric(i), "start_break_pos"] # select the first coefficient 
    end_break   <- start_break + floor(mot_len[mot_len[,1]==as.numeric(i),2] / dist_knots) + mot_order - 2 
    for (n in seq_along(start_break)) {
      base_curve_coeff[start_break[n]:end_break[n]] <- mot_coeff[[i]]
      if(!only_der) {
        base_curve_coeff[start_break[n]:end_break[n]] <- base_curve_coeff[start_break[n]:end_break[n]] + rep(runif(1,min=coeff_min_shift,max=coeff_max_shift),length(mot_coeff[[i]])) 
      }
    }
  }
  
  # Errorless curve
  # create fd object
  # Convert base_curve_coeff and basis to fd object
  fd_curve <- fd(base_curve_coeff, base_curve$basis)
  
  # Generate curve vector with no noise
  or_y <- generate_curve_vector(fd_curve = fd_curve)
  
  background <- list("or_coeff" = base_curve$or_coeff,
                     "no_error_y" = base_y)
  no_error_res <- list("or_coeff" = base_curve_coeff,
                       "motif_y" = or_y)
  # Error curve
  # add extra error on motif
  err_y_mat <- list()
  noise_str <- lapply(noise_str,.transform_to_matrix)
  SNR <- vector("list",nrow(noise_str[[1]]))
  for(k in 1:nrow(noise_str[[1]])) {
    err_y <- or_y
    SNR_num <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
    SNR_den <- data.frame(xmin = numeric(), xmax = numeric(), SNR = numeric(), stringsAsFactors = FALSE)
    for(i in names(mot_weights)){
      error_str_i_k <- noise_str[[as.numeric(i)]][k,]
      start_break <- mot_pattern[mot_pattern[,1]==as.numeric(i), "start_break_pos"] # select the first coefficient 
      start_point <- (start_break-1)*dist_knots
      end_point   <- start_point + mot_len[mot_len[,1]==as.numeric(i),"len"]
      err_y <- add_error_to_motif(err_y, error_str_i_k, start_point, end_point,k)
      for(n in seq_along(start_point)) {
        SNR_num <- rbind(SNR_num,data.frame(xmin =start_point[n],xmax = end_point[n],SNR = var(no_error_res$motif_y[start_point[n]:end_point[n]])))
      }
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
    for(i in names(mot_weights)) {
      start_break <- mot_pattern[mot_pattern[,1]==as.numeric(i), "start_break_pos"] # select the first coefficient 
      start_point <- (start_break-1)*dist_knots
      end_point   <- start_point + mot_len[mot_len[,1]==as.numeric(i),"len"]
      for(n in seq_along(start_point)) {
        SNR_den <- rbind(SNR_den,data.frame(xmin=start_point[n],xmax=end_point[n],SNR=var(yy[start_point[n]:end_point[n]] - no_error_res$motif_y[start_point[n]:end_point[n]])))
      }
    }
    SNR[[k]] <- SNR_num
    SNR[[k]]$SNR <- 10 * log10(SNR_num$SNR / SNR_den$SNR) #transform SNR in decibel
    err_y_mat[[k]] <- yy
  }

  error_res <- list("noise_structure" = noise_str,
                    "noise_y" = err_y_mat) 
  res <- list("basis" = base_curve$basis,
              "background" = background,
              "no_noise" = no_error_res,
              "with_noise" = error_res,
              "SNR" = SNR) 
  # create list with 
  return(res)
}

.generate_coefficients <- function(motif_i, distrib, dist_knots, norder, coeff_min, coeff_max) {
  # Calculate the length of the coefficients vector
  l <- motif_i$len / dist_knots + norder - 1
  if(is.numeric(distrib)) {
    motif_i$coeffs <- sample(distrib,size = l,replace = TRUE)
  }
  else if (distrib == "unif") {
    # Generate coefficients from a uniform distribution
    motif_i$coeffs <- runif(l, min = coeff_min, max = coeff_max)
  } else if (distrib == "beta") {
    # Generate coefficients from a beta distribution and scale to the desired range
    motif_i$coeffs <- coeff_min + rbeta(l, shape1 = 0.45, shape2 = 0.45) * (coeff_max - coeff_min)
  } else {
    stop("Wrong 'distrib': ", distrib)
  }
  return(motif_i)
}

# Check for overlaps and fitting within each curve
.check_fits <- function(df,min_dist_motifs,len) {
  df <- df[order(df$start), ]
  # Check for overlaps
  for (i in seq_len(nrow(df) - 1)) {
    if (df$end[i] + min_dist_motifs > df$start[i + 1]) {
      return(FALSE)
    }
  }
  
  # Check if the last motif fits within the curve length
  if (nrow(df) > 0 && df$end[nrow(df)] > len) {
    return(FALSE)
  }
  
  return(TRUE)
}

.resample <- function(x, ...) x[sample.int(length(x), ...)]



