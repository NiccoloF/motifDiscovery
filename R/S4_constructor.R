#' @title Create motifSimulation Object
#' @description It represents the constructor of the S4 class 'motifSimulation'..
#' @param N An integer specifying the number of background curves to be generated (mandatory).
#' @param len An integer specifying the length of the background curves (mandatory).
#' @param mot_details A list outlining the definitions of the motifs to be included. Each motif is characterized by its length, a set of coefficients that may be optionally specified, and the number of occurrences. These occurrences can be indicated either by specific positions within the curves or by a total count. In the latter case, the algorithm will randomly position the motifs throughout the curves (mandatory).
#' @param norder An integer specifying the order of the B-splines (default = 3).
#' @param coeff_min Additive coefficients to be incorporated into the generation of coefficients for the background curves (default = -15).
#' @param coeff_max Additive coefficients to be incorporated into the generation of coefficients for the background curves (default = 15).
#' @param dist_knots An integer specifying the distance between two consecutive knots (default = 10).
#' @param min_dist_motifs An integer specifying the minimum distance between two consecutive motifs embedded in the same curve (default = 'norder' * 'dist_knots').
#' @param distribution A character string specifying the distribution from which the coefficients of the background curves are generated. You can choose between a uniform distribution or a beta distribution. Alternatively, you can pass a vector representing the empirical distribution from which you wish to sample (default = "unif").
#' @return An object of class motifSimulation
#' @export
#' @examples
#' \donttest{
#' mot_len <- 100
#' motif_str <- rbind.data.frame(c(1, 1, 20),
#'                               c(2, 1, 2), 
#'                               c(1, 3, 1),
#'                               c(1, 2, 1),
#'                               c(1, 2, 15),
#'                               c(1, 4, 1),
#'                               c(2, 5, 1),
#'                               c(2, 7, 1),
#'                               c(2,17,1))
#' 
#' names(motif_str) <- c("motif_id", "curve","start_break_pos")
#' 
#' mot1 <- list("len" = mot_len, #length
#'              "coeffs" = NULL, # weights for the motif
#'              "occurrences" = motif_str %>% filter(motif_id == 1))
#' 
#' mot2 <- list("len" = mot_len,
#'              "coeffs" = NULL,
#'              "occurrences" = motif_str %>% filter(motif_id == 2))
#' 
#' mot_details <- list(mot1,mot2)
#' 
#' # MATRIX ERROR 
#' noise_str <- list(rbind(rep(2, 100)),
#'                   rbind(rep(0.0, 100)))
#' 
#' builder <- funMoDisco::motifSimulationBuilder(N = 20,len = 300,mot_details,
#'                                         distribution = 'beta')
#' }

motifSimulationBuilder <- function(N,len,mot_details,norder = 3,
                                   coeff_min=-15,coeff_max=15,
                                   dist_knots=10,min_dist_motifs=NULL,
                                   distribution = 'unif') {
  ########################## CHECKS #################################
  if(is.null(min_dist_motifs)) {
    min_dist_motifs <- norder * dist_knots
  }
  # check N, dist_knots and len
  if((N%%1!=0)|(N<1))
    stop('Invalid \'N\'.')
  if(length(mot_details) == 0) {
    return(new("motifSimulation",N = N,
               mot_details = list(),
               motifs_in_curves = list(),
               distribution = distribution,
               dist_knots=dist_knots,len=len,norder=norder,
               coeff_min=coeff_min,coeff_max=coeff_max,
               min_dist_motifs=min_dist_motifs))
  }
  if((dist_knots%%1!=0)|(dist_knots<1))
    stop('Invalid \'dist_knots\'.')
  if(TRUE %in% ((len%%1!=0)|(len<1)|(sum(len%%dist_knots)>0)))
    stop('Invalid \'len\'.')
  mot_details <- lapply(mot_details,function(mot){names(mot) <- c("len","coeffs","occurrences")
                                                  if(is.data.frame(mot[[3]]) || is.matrix(mot[[3]])) { 
                                                    names(mot[[3]]) <- c("motif_id", "curve","start_break_pos") }
                                                  return(mot)})

  weights_defined <- sapply(mot_details, function(x) !is.null(x$coeffs))
  if (!(all(weights_defined) || all(!weights_defined))) {
    stop("Inconsistent coefficients field: some elements have it defined, others do not.")
  }
  
  len_motifs <- sapply(mot_details,function(x){x$len},simplify = TRUE)
  nmotifs <- length(mot_details)
  if(TRUE %in% weights_defined) {
    for(i in 1:nmotifs){
      nbasis <- len_motifs[i]/dist_knots + norder - 1
      mot_details[[i]]$coeffs <- mot_details[[i]]$coeffs[1:nbasis]
      
    }
  }
  # check nmotifs and len_motifs
  if((nmotifs%%1!=0)|(nmotifs<1))
    stop('Invalid \'nmotifs\'.')
  if(TRUE %in% ((len_motifs%%1!=0)|(len_motifs<1)|(sum(len_motifs%%dist_knots)>0)|(TRUE %in% (len_motifs>len))))
    stop('Invalid \'len_motifs\'.')
  
  if (FALSE %in% unlist(lapply(mot_details, function(x) {
    if (all(weights_defined)) {
      x$len == (length(x$coeffs)-norder+1)*dist_knots
    } else {
      TRUE
    }
  }))) {
    stop("The length of the motifs is incompatible")
  }
  # check min_dist_motifs
  if((min_dist_motifs < norder*dist_knots)|(min_dist_motifs%%dist_knots>0))
    stop('Invalid \'min_dist_motifs\'.')
  # check freq_motifs and convert it into a list with the frequencies for the different motifs in the different curves
  
  freq_motifs_vector <- sapply(mot_details,function(mot_i){ifelse(is.data.frame(mot_i$occurrences),dim(mot_i$occurrences['curve'])[1],NA)})
  is_appearance_defined <- !(NA %in% freq_motifs_vector)
  if(sum(is.na(freq_motifs_vector)) != 0 && sum(is.na(freq_motifs_vector)) != nmotifs) {
    stop("Inconsistent appearance field: some elements have it defined, others do not.")
  }
  if(is_appearance_defined) {
    freq_motifs_vector <- freq_motifs_vector[!is.na(freq_motifs_vector)]
    if(TRUE %in% ((freq_motifs_vector%%1!=0)|(freq_motifs_vector<1)|(sum(norder-1+rep(len/dist_knots,length.out=N)-2*(norder-1))<(sum(rep(len_motifs/dist_knots+norder-1,length.out=nmotifs)*rep(freq_motifs_vector,length.out=nmotifs))+max(0,sum(rep(freq_motifs_vector, length.out = nmotifs))-N)*(min_dist_motifs/dist_knots-norder+1)))))
      stop('Invalid \'freq_motifs\'.')
    motif_str_list <- lapply(mot_details, function(x) {
      df <- x$occurrences
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
                            .check_fits(sub_data, min_dist_motifs = min_dist_motifs,len = len)
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
        curves <- as.vector(sample(N,mot_details[[motif_j]]$occurrences,replace=TRUE))
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
        count <- sum(mot_details[[motif_j]]$occurrences$curve == curve_i)
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
      mot_details[[mot]]$occurrences <- data.frame(motif_id = integer(0),
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
        mot_details[[mot]]$occurrences <<- rbind(mot_details[[mot]]$occurrences, df)
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
             min_dist_motifs=min_dist_motifs))
}
