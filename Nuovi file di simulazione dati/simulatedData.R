library(fda)
library(tidyverse)
library(ggplot2)
library(scales)

#
#mot_details <- list("mot_len" = mot_len, #length
#                    "mot_order" = mot_order, # order of the bspline
#                    "weight" = rbeta(300, 1, 3), # weights for the motif
#                    "appearance" = motif_str %>% filter(motif_id == 1)) #pattern
#                    "distr" = NULL

#motif_str <- rbind.data.frame(c(1, 1, 10),
#                              c(1, 1, 35),
#                              c(1, 1, 89),
#                              c(1, 2, 15),
#                              c(1, 2, 30))

norder <- 3
mot_len <- 100
mot_order <-  norder
mot_weights <- rbeta(300, 1, 3) #runif(300, -15, 15)

# Generate multiple curves with a motif structure imposed
N <- 20

#motif_str <- rbind.data.frame(c(1, 1, 10),
#                              c(1, 1, 35),
#                              c(1, 1, 89),
#                              c(1, 2, 15),
#                              c(1, 2, 30),
#                              c(1, 4, 5),
#                              c(1, 4, 60),
#                              c(1, 5, 89),
#                              c(1, 7, 70),
#                              c(1,10, 18),
#                              c(1,10, 32),
#                              c(1,12, 40),
#                              c(1,12, 21),
#                              c(1,15, 17),
#                              c(1,17, 81),
#                              c(1,17, 20),
#                              c(1,18, 3),
#                              c(1,20, 20),
#                              c(1,20, 69),
#                              c(2,17,1),
#                              c(2,17,55))

motif_str <- rbind.data.frame(c(1, 1, 20),
                              c(1, 1, 2), #9 
                              c(1, 3, 1),
                              c(1, 2, 1),
                              c(1, 2, 15),
                              c(1, 4, 1),
                              c(1, 5, 1),
                              c(1, 7, 1),
                              c(2,17,1))

names(motif_str) <- c("motif_id", "curve","start_break_pos")
dim(motif_str)

mot1 <- list("len" = mot_len, #length
             "weights" = NULL, # weights for the motif
             "appearance" = 5) #pattern motif_str %>% filter(motif_id == 1)

mot2 <- list("len" = mot_len,
             "weights" = NULL,
             "appearance" = 6)

mot_details <- list(mot1,mot2)

setClass("motifSimulation",
         slots = list(
           N = "numeric",            # Number of curves
           mot_details = "list",     # Details of motifs 
           motifs_in_curves = "list",# Details of curves
           distribution = "character",  # distribution of the weights
           dist_knots = "numeric",   # Distance between knots 
           len = "numeric",          # Length of the curve or curves
           norder = "numeric",       # Order of spline 
           coeff_min = "numeric",    # Minimum coefficient value
           coeff_max = "numeric",    # Maximum coefficient value
           min_dist_motifs = "numeric", # Minimum distance between motifs
           is_appearance_defined = "logical" # boolean
         ))


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
                   "but the motif is long",length(start_point:end_point),". -- Extension"))
    # Last value
    last_value <- tail(error_str, 1)
    
    # Extention 
    error_str <- c(error_str, rep(last_value, length(start_point:end_point) - length(error_str)))
  }else if(length(start_point:end_point) - 1 < length(error_str)) {
    
    warning(paste("error structure on row",k,"is long",length(error_str),
                   "but the motif is long",length(start_point:end_point),". -- Truncation"))
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
  
  # # Create a data frame for plotting
  # df <- data.frame(x = 1:length(or_y), 
  #                  or_y = or_y,
  #                  no_error_y = base_curve$no_error_y)
  # colnames(df) <- c("x","or_y","no_error_y")
  
  # Plot using ggplot2
  # p <- ggplot(df, aes(x = x)) +
  #   geom_line(aes(y = no_error_y, color = "Original Curve"), linetype = "solid") +
  #   geom_line(aes(y = or_y, color = "Curve with Motif"), linetype = "solid") +
  #   labs(title = 'Original vs Curve with Motif without noise', x = 'X-axis', y = 'Y-axis') +
  #   theme_minimal() +
  #   scale_color_manual(values = c("Original Curve" = "black", "Curve with Motif" = "grey60")) +
  #   guides(color = guide_legend(override.aes = list(linetype = "solid"))) +
  #   theme(legend.position = "top") +
  #   theme(legend.key = element_blank(),  # Remove legend key (box)
  #         legend.title = element_blank(),  # Remove legend title
  #         legend.text = element_text(size = 12)) +  # Adjust legend text size
  #   labs(color = NULL)  # Remove color legend title
  # no_error_res <- list("or_coeff" = base_curve_coeff,
  #                      "basis" = base_curve$basis,
  #                      "no_error_y" = or_y, # curve as vector with no error
  #                      "sd_noise_level" = 0,
  #                      "motif_info" = motif_info)
  # 
  # print(p)
  
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
      err_y <- add_error_to_motif(err_y, error_str_k, start_point, end_point,k)
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
  #
  # #  plot curves without noise
  # x <- seq_along(or_y)
  # 
  # # Create a data frame for ggplot2
  # df <- data.frame(x = x, or_y = or_y)
  # 
  # # Plot using ggplot2
  # p <- ggplot(df, aes(x = x, y = or_y)) +
  #   geom_line(aes(color = "Curva"), linewidth = 1) +
  #   labs(title = "Curves with Motifs Embedded", x = "X-axis", y = "Y-axis") +
  #   theme_minimal() +
  #   theme(legend.position = "top") +  # Posizione della legenda in alto
  #   scale_color_manual(values = c("Curva" = "black", "Motif" = "red"), 
  #                      labels = c("Curva", "Motif"))  # Personalizzazione dei colori e delle etichette
  # 
  # # Add motifs embedded
  # for (i in 1:nrow(motif_recap)) {
  #   motif_df <- data.frame(x = x[motif_recap[i, "start_point"]:motif_recap[i, "end_point"]],
  #                          y = or_y[motif_recap[i, "start_point"]:motif_recap[i, "end_point"]])
  #   
  #   p <- p + geom_line(data = motif_df, aes(x = x, y = y, color = "Motif"), linewidth = 1.5)
  # }
  # # Display the plot
  # print(p)
  
  error_res <- list("error_structure" = error_str,
                    "error_y" = err_y_mat) 
  res <- list("no_error" = no_error_res,
              "with_error" = error_res)
  # create list with 
  return(res)
}



motifSimulationBuilder <- function(curve_details = list(N=20,
                                                        len = 300,
                                                        norder = 3,
                                                        coeff_min=-15,
                                                        coeff_max=15,
                                                        dist_knots=10,
                                                        min_dist_motifs=norder * dist_knots),
                                   mot_details,
                                   distribution = 'unif') {
  ########################## UNPACKING #################################
  N <- NULL
  len <- NULL
  norder <- NULL
  coeff_min <- NULL
  coeff_max <- NULL
  dist_knots <- NULL
  min_dist_motifs <- NULL

  if('N' %in% names(curve_details)) {
    N <- curve_details[['N']]
  }else{
    N <- curve_details[[1]]
  }
  if('len' %in% names(curve_details)) {
    len <- curve_details[['len']]
  }else{
    len <- curve_details[[2]]
  }
  if('norder' %in% names(curve_details)) {
    norder <- curve_details[['norder']]
  }else{
    norder <- curve_details[[3]]
  }
  if('coeff_min' %in% names(curve_details)) {
    coeff_min <- curve_details[['coeff_min']]
  }else{
    coeff_min <- curve_details[[4]]
  }
  if('coeff_max' %in% names(curve_details)) {
    coeff_max <- curve_details[['coeff_max']]
  }else{
    coeff_max <- curve_details[[5]]
  }
  if('dist_knots' %in% names(curve_details)) {
    dist_knots <- curve_details[['dist_knots']]
  }else{
    dist_knots <- curve_details[[6]]
  }
  if('min_dist_motifs' %in% names(curve_details)) {
    min_dist_motifs <- curve_details[['min_dist_motifs']]
  }else{
    min_dist_motifs <- curve_details[[7]]
  }

  ########################## AUXILIARY FUNCTIONS #################################
  # Function to generate coefficients based on the distribution
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
  .check_fits <- function(df) {
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
  
  ########################## CHECKS #################################
  # check N, dist_knots and len
  if((N%%1!=0)|(N<1))
    stop('Invalid \'N\'.')
  if((dist_knots%%1!=0)|(dist_knots<1))
    stop('Invalid \'dist_knots\'.')
  if(TRUE %in% ((len%%1!=0)|(len<1)|(sum(len%%dist_knots)>0)))
    stop('Invalid \'len\'.')
  if(length(mot_details) == 0)
    stop('Invalid \'mot_details\'.')
  
  weights_defined <- sapply(mot_details, function(x) !is.null(x$weights))
  if (!(all(weights_defined) || all(!weights_defined))) {
    stop("Inconsistent weights field: some elements have it defined, others do not.")
  }
  
  len_motifs <- sapply(mot_details,function(x){x$len},simplify = TRUE)
  nmotifs <- length(mot_details)
  if(TRUE %in% weights_defined) {
    for(i in 1:nmotifs){
      nbasis <- len_motifs[i]/dist_knots + norder - 1
      mot_details[[i]]$weights <- mot_details[[i]]$weights[1:nbasis]
      
    }
  }
  # check nmotifs and len_motifs
  if((nmotifs%%1!=0)|(nmotifs<1))
    stop('Invalid \'nmotifs\'.')
  if(TRUE %in% ((len_motifs%%1!=0)|(len_motifs<1)|(sum(len_motifs%%dist_knots)>0)|(TRUE %in% (len_motifs>len))))
    stop('Invalid \'len_motifs\'.')
  
  if (FALSE %in% unlist(lapply(mot_details, function(x) {
    if (all(weights_defined)) {
      x$len == (length(x$weights)-norder+1)*dist_knots
    } else {
      TRUE
    }
  }))) {
    stop("The length of the motifs is incompatible ")
  }
  # check min_dist_motifs
  if((min_dist_motifs < norder*dist_knots)|(min_dist_motifs%%dist_knots>0))
    stop('Invalid \'min_dist_motifs\'.')
  # check freq_motifs and convert it into a list with the frequencies for the different motifs in the different curves
  
  freq_motifs_vector <- sapply(mot_details,function(mot_i){ifelse(is.data.frame(mot_i$appearance),dim(mot_i$appearance['curve'])[1],NA)})
  is_appearance_defined <- !(NA %in% freq_motifs_vector)
  if(sum(is.na(freq_motifs_vector)) != 0 && sum(is.na(freq_motifs_vector)) != nmotifs) {
    stop("Inconsistent appearance field: some elements have it defined, others do not.")
  }
  if(is_appearance_defined) {
    freq_motifs_vector <- freq_motifs_vector[!is.na(freq_motifs_vector)]
    if(TRUE %in% ((freq_motifs_vector%%1!=0)|(freq_motifs_vector<1)|(sum(norder-1+rep(len/dist_knots,length.out=N)-2*(norder-1))<(sum(rep(len_motifs/dist_knots+norder-1,length.out=nmotifs)*rep(freq_motifs_vector,length.out=nmotifs))+max(0,sum(rep(freq_motifs_vector, length.out = nmotifs))-N)*(min_dist_motifs/dist_knots-norder+1)))))
      stop('Invalid \'freq_motifs\'.')
    motif_str_list <- lapply(mot_details, function(x) {
      df <- x$appearance
      len <- x$len
      df$len <- len  # Add the length to each motif dataframe
      return(df)
    })
    motif_str_new <- do.call(rbind, motif_str_list)
    
    # Convert columns to numeric if necessary
    motif_str_new[] <- lapply(motif_str_new, as.numeric)
    
    # Add columns for start and end positions
    motif_str_new <- motif_str_new %>%
      mutate(
        start = (start_break_pos - 1) * dist_knots,
        end = start + len
      )
    # Apply the check for each curve
    valid_curves <- sapply(split(motif_str_new, motif_str_new$curve), .check_fits)
    
    # Check if all curves have valid motifs
    if(FALSE %in% valid_curves) {
      stop('Invalid position of the knots.')
    }
  }
  
  
  freq_check=1
  it=0
  
  # Compute the frequency(if not defined) of the motifs inside the curves 
  if(!is_appearance_defined) {
    while((sum(freq_check)>0)&&(it<10000)){
      it=it+1
      freq_motifs <- matrix(data = 0,nrow = N,ncol = nmotifs)
      for(motif_j in 1:nmotifs) {
          curves <- as.vector(sample(N,mot_details[[motif_j]]$appearance,replace=TRUE))
          for (curve in curves) {
            freq_motifs[curve, motif_j] <- freq_motifs[curve, motif_j] + 1
            }
        }
        freq_motifs=split(freq_motifs,rep(1:N,nmotifs))# construct N list(for each curve). Each list_j contains how many times the motif_i appears in the curve j
        freq_check=mapply(function(freq_motifs_i,len_i) ((len_i-2*(norder-1))<(sum(rep(len_motifs/dist_knots+norder-1,length.out=nmotifs)*rep(freq_motifs_i,length.out=nmotifs))+(sum(freq_motifs_i)-1)*(min_dist_motifs/dist_knots-norder+1))),
                        freq_motifs,norder-1+rep(len/dist_knots,length.out=N))
    }
  if(it==10000)
    stop('Tried 10.000 random configurations of motif frequencies in the different curves, unable to find a valid configuration. Please select lower \'freq_motifs\'.')
  }else {
    freq_motifs <- matrix(data = 0,nrow = N,ncol = nmotifs)
    for(motif_j in 1:nmotifs) {
      for(curve_i in 1:N){
        count <- sum(mot_details[[motif_j]]$appearance$curve == curve_i)
        freq_motifs[curve_i, motif_j] <-  freq_motifs[curve_i, motif_j] + count
      }
    }
    ifelse(TRUE %in% mapply(function(freq_motifs_i,len_i) ((len_i-2*(norder-1))<(sum(rep(len_motifs/dist_knots+norder-1,length.out=nmotifs)*rep(freq_motifs_i,length.out=nmotifs))+(sum(freq_motifs_i)-1)*(min_dist_motifs/dist_knots-norder+1))),
                            freq_motifs,norder-1+rep(len/dist_knots,length.out=N)),"",stop('Please select lower \'freq_motifs\'.'))
  }
  
  ########################## POSITION ASSIGNMENT ###############################
  # randomly assign motif positions in curves if they are not already assigned
  # At this point we only have the coefficients of the motifs
  motifs_in_curves <- vector("list", N)
  if (!is_appearance_defined) {
    # Initialize appearance for each motif
    for (mot in 1:nmotifs) {
      mot_details[[mot]]$appearance <- data.frame(motif_id = integer(0),
                                                  curve = integer(0),
                                                  start_break_pos = integer(0),
                                                  coeff_pos = integer(0))
    }
    
    motifs_in_curves <- mapply(function(freq_motifs_i, len_i, curve_i) {
      if (sum(freq_motifs_i) == 0) {
        return(NULL) # No motifs embedded in this curve
      }
      
      id_motifs <- .resample(rep(seq_along(freq_motifs_i), freq_motifs_i))
      len_elements <- c(norder - 1, rep(len_motifs / dist_knots + norder - 1, length.out = nmotifs)[id_motifs] + c(rep((min_dist_motifs / dist_knots - norder + 1), sum(freq_motifs_i) - 1), 0), norder - 1)
      gaps_tot <- len_i - sum(len_elements) # number of free coefficients
      gaps <- diff(c(0, sort(sample(gaps_tot + sum(freq_motifs_i), sum(freq_motifs_i))))) - 1
      coeff_pos <- cumsum(len_elements[seq_along(id_motifs)]) + 1 + cumsum(gaps) # Starting coefficient
      
      for (mot in unique(id_motifs)) {
        indices <- which(id_motifs == mot)
        df <- data.frame(
          motif_id = mot,
          curve = rep(curve_i, length(indices)),
          start_break_pos = coeff_pos[indices] # coincide with the coeffs position
        )
        mot_details[[mot]]$appearance <<- rbind(mot_details[[mot]]$appearance, df)
      }
      
      return(list(motif_id=id_motifs,starting_coeff_pos=coeff_pos)) # returns N list(one for each curve) with the id of the motif embedded and the starting position
    }, freq_motifs, norder - 1 + rep(len / dist_knots, length.out = N), 1:N,SIMPLIFY = FALSE)
  } else {
    splited_curves <- split(motif_str_new, motif_str_new$curve)
    names(motifs_in_curves) <- as.character(1:N)
    for (name in as.character(1:N)) {
      if (name %in% names(splited_curves)) {
        motifs_in_curves[[name]] <- list(motif_id = splited_curves[[name]]$motif_id,starting_coeff_pos = splited_curves[[name]]$start_break_pos)
      }
    }
  }
  
  ########################## WEIGHTS GENERATION #################################
  # if weights are not provided, randomly generate them according to the distribution
  if(all(!weights_defined)) {
    mot_details <- mapply(.generate_coefficients,mot_details,
                           MoreArgs = list(distrib = distribution,
                                           dist_knots = dist_knots,
                                           norder = norder, 
                                           coeff_min = coeff_min,
                                           coeff_max = coeff_max),
                          SIMPLIFY = FALSE)
  }
  
  ########################## RETURN #################################
  return(new("motifSimulation",N = N,
             mot_details = mot_details,
             motifs_in_curves = motifs_in_curves,
             distribution = distribution,
             dist_knots=dist_knots,len=len,norder=norder,
             coeff_min=coeff_min,coeff_max=coeff_max,
             min_dist_motifs=min_dist_motifs,
             is_appearance_defined = is_appearance_defined))
}

curve_details = list(N=20,
                     len = 300,
                     norder = 3,
                     coeff_min=-15,
                     coeff_max=15,
                     dist_knots=10,
                     min_dist_motifs=30)
b <- motifSimulationBuilder(curve_details,
                            mot_details,
                            distribution = 'unif')



# error_str can be a vector or a matrix with n_cols = motif_len
setGeneric("generateCurves", function(object,error_str) {
  standardGeneric("generateCurves")
})

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
          or_y_no_error <- generate_curve_vector(fda_no_error)
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
        or_y_no_error <- generate_curve_vector(fda_no_error)
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
        or_y <- Map(generate_curve_vector,fda_with_error)
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
      generate_background_curve(object@len, object@dist_knots, object@norder,coeff, add_noise = TRUE)
    })
    for(i in 1:length(object@mot_details)){
      print(paste("Dealing with motif", i))
      temp_mot <- object@mot_details[[i]]
      curve_ids <- unique(temp_mot$appearance$curve)
      for(j in curve_ids){
        print(paste(" --- Adding motif", i, "to curve", j))
        temp_curve <- fd_curves[[j]]
        temp_pattern <- (temp_mot$appearance %>% filter(curve == j)) %>% dplyr::select(motif_id, start_break_pos)
        fd_curves[[j]] <- add_motif(base_curve  = temp_curve,
                                    mot_pattern = temp_pattern,
                                    mot_len     = temp_mot$len,
                                    dist_knots  = object@dist_knots,
                                    mot_order   = object@norder,
                                    mot_weights = temp_mot$weight,
                                    error_str   = error_str)
      }
    }
    .transform_list <- function(lst) {
      # Looping the list
      for (i in seq_along(lst)) {
        # check"with_error" field
        if (!"with_error" %in% names(lst[[i]])) {
          # Create a new sub-list
          lst[[i]] <- list(
            no_error = lst[[i]],
            with_error = list(
              error_structure = .transform_to_matrix(error_str),
              error_y = NULL)
          )
        }
      }
      return(lst)
    }
    fd_curves <- .transform_list(fd_curves)
  }
  return(fd_curves = fd_curves)
})

#############   
error_str <- rbind(
  rep(2, 100), # constant and identical
  c(rep(0.1,50), rep(2, 50)), # sd 0.1 first, sd 1 later
  c(rep(2, 50), rep(0.1, 50)), # sd 1 first, 0.1 later
  c(seq(2, 0.1, len = 50), rep(0.1, 50)) # decreasing first, constant later
)
##############  
curves <- generateCurves(b,c(0.1,0.5))


setGeneric("plot", function(object,curves,path) 
  standardGeneric("plot")
)

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


plot(b,curves,"plot_simulated_data")

