#' @title Create MyS4Class Object
#' @description A constructor function for the S4 class MyS4Class.
#' @param name A character string representing the name.
#' @param value A numeric value.
#' @return An object of class MyS4Class.
#' @export
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
    valid_curves <- mapply(function(sub_data) {
                            mtfd:::.check_fits(sub_data, min_dist_motifs = min_dist_motifs,len = len)
                          },
                          split(motif_str_new, motif_str_new$curve), 
                          SIMPLIFY = TRUE)
    
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
    freq_motifs=split(freq_motifs,rep(1:N,nmotifs))
    ifelse(TRUE %in% mapply(function(freq_motifs_i,len_i) ((len_i-2*(norder-1))<(sum(rep(len_motifs/dist_knots+norder-1,length.out=nmotifs)*rep(freq_motifs_i,length.out=nmotifs))+(sum(freq_motifs_i)-1)*(min_dist_motifs/dist_knots-norder+1))),
                            freq_motifs,norder-1+rep(len/dist_knots,length.out=N)),stop('Please select lower \'freq_motifs\'.'),"")
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
      
      id_motifs <- mtfd:::.resample(rep(seq_along(freq_motifs_i), freq_motifs_i))
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
    mot_details <- mapply(mtfd:::.generate_coefficients,mot_details,
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