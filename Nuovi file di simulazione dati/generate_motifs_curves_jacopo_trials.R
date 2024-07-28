library(fda)
library(tidyverse)
library(ggplot2)

generate_curve_vector <- function(fd_curve, step_by = 1, Lfdobj = 0){
  # from an fd object generate the curve as a vector of len fixed by the fd_obj
  # with step equal to step_by. It generates the derivative Lfdobj
  x <- seq(0, fd_curve$basis$rangeval[2], by = step_by)
  y <- eval.fd(x, fd_curve, Lfdobj = Lfdobj)  # return the evaluation of the Lfdobj derivative
  return(y)
}

# add an additive noise to the final curve(motif)
add_error_to_motif <- function(or_y, error_str, start_point, end_point){
  
  # --- additive noise
  # noise to be added as a percentage of sd(motif) with mean 0
  noise_coeff <- rnorm(length(start_point:end_point),
                       mean = 0,
                       sd = sd(or_y[start_point:end_point])*error_str)
  # add noise
  err_y <- or_y
  err_y[start_point:end_point] <- err_y[start_point:end_point] + noise_coeff
  
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
              "no_error_y" = or_y, # curve as vector with no error
              "sd_noise_level" = 0,
              "motif_info" = NULL) # level of error used
  return(res)
}




len <- 1010
dist_knots <- 10
norder <- 3
weights <- rbeta(300, 1, 3)
base_curve   <- generate_background_curve(len, dist_knots, norder, weights, add_noise = TRUE)
plot(base_curve$no_error_y, type='l')

# add motif to background curve (base_curve) passing many info about the motif
# such as the pattern, the length, the distance of the knots in the bspline,
# the motif order, the motif weight
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
  motif_info <- base_curve$motif_info
  base_curve_coeff <- base_curve$or_coeff
  motif_recap <- cbind.data.frame()
  err_y_list <- list()
  add_err_y_list <- list()
  for(i in 1:nrow(mot_pattern)){
    start_break <- mot_pattern[i, "start_break"] # select the first coefficient 
    start_point <- start_break*dist_knots # select the first knots
    end_break   <- start_break + mot_basis$nbasis - 1 # select the last coefficient 
    end_point   <- end_break*dist_knots # select the last knots
    motif_recap[i, "motif_id"]   <- mot_pattern[i, "motif_id"]
    motif_recap[i,"start_break"] <- start_break
    motif_recap[i,"start_point"] <- start_point
    motif_recap[i,"end_break"]   <- end_break
    motif_recap[i,"end_point"]   <- end_point
    motif_info <- rbind.data.frame(motif_info, motif_recap) %>% distinct() #(!!!) find a way to avoid distinct(). Keeo unique rows
    base_curve_coeff[start_break:end_break] <- mot_coeff
  }
  
  # Errorless curve
  # create fd object
  # Convert base_curve_coeff and basis to fd object
  fd_curve <- fd(base_curve_coeff, base_curve$basis)
  
  # Generate curve vector with no noise
  or_y <- generate_curve_vector(fd_curve = fd_curve)
  
  # Create a data frame for plotting
  df <- data.frame(x = 1:length(or_y), 
                   or_y = or_y,
                   no_error_y = base_curve$no_error_y)
  colnames(df) <- c("x","or_y","no_error_y")

  # Plot using ggplot2
  p <- ggplot(df, aes(x = x)) +
    geom_line(aes(y = no_error_y, color = "Original Curve"), linetype = "solid") +
    geom_line(aes(y = or_y, color = "Curve with Motif"), linetype = "solid") +
    labs(title = 'Original vs Curve with Motif without noise', x = 'X-axis', y = 'Y-axis') +
    theme_minimal() +
    scale_color_manual(values = c("Original Curve" = "black", "Curve with Motif" = "grey60")) +
    guides(color = guide_legend(override.aes = list(linetype = "solid"))) +
    theme(legend.position = "top") +
    theme(legend.key = element_blank(),  # Remove legend key (box)
          legend.title = element_blank(),  # Remove legend title
          legend.text = element_text(size = 12)) +  # Adjust legend text size
    labs(color = NULL)  # Remove color legend title
   no_error_res <- list("or_coeff" = base_curve_coeff,
                       "basis" = base_curve$basis,
                       "no_error_y" = or_y, # curve as vector with no error
                       "sd_noise_level" = 0,
                       "motif_info" = motif_info) 
   
   print(p)
  
  # Error curve
  # add extra error on motif
  
  if(!is.null(error_str)){
    err_y_mat <- list()
    for(k in 1:nrow(error_str)){
      err_y <- or_y
      error_str_k <- error_str[k,]
      for(i in 1:nrow(mot_pattern)){
        start_break <- mot_pattern[i, "start_break"]
        start_point <- start_break*dist_knots
        end_break   <- start_break + mot_basis$nbasis - 1 # mot_basis$nbasis - order + 2 
        end_point   <- end_break*dist_knots
        err_y <- add_error_to_motif(err_y, error_str_k, start_point, end_point)
      }

      # and then smooth it again
      mot_breaks <- seq(from = 1, length(err_y), by = dist_knots)
      basisobj = create.bspline.basis(norder = 3, breaks = mot_breaks)
      #plot(basisobj)
      ys = smooth.basis(argvals = 1:(length(err_y)),
                        y = err_y,
                        fdParobj = basisobj) # smooth again given the error on data(pointwise)
      # smoothed again after adding error
      yy <- eval.fd(1:(length(err_y)), ys$fd)
      err_y_mat[[k]] <- yy
    }
      #lines(yy, col='blue')
    
  error_res <- list("error_structure" = error_str,
                    "error_y" = err_y_mat)
  
  res <- list("no_error" = no_error_res,
              "with_error" = error_res)
  }
#
#  plot curves without noise
   x <- seq_along(or_y)
   
   # Create a data frame for ggplot2
   df <- data.frame(x = x, or_y = or_y)
   
   # Plot using ggplot2
   p <- ggplot(df, aes(x = x, y = or_y)) +
     geom_line(aes(color = "Curva"), size = 1) +
     labs(title = "Curves with Motifs Embedded", x = "X-axis", y = "Y-axis") +
     theme_minimal() +
     theme(legend.position = "top") +  # Posizione della legenda in alto
     scale_color_manual(values = c("Curva" = "black", "Motif" = "red"), 
                        labels = c("Curva", "Motif"))  # Personalizzazione dei colori e delle etichette
   
   # Add motifs embedded
   for (i in 1:nrow(motif_recap)) {
     motif_df <- data.frame(x = x[motif_recap[i, "start_point"]:motif_recap[i, "end_point"]],
                            y = or_y[motif_recap[i, "start_point"]:motif_recap[i, "end_point"]])
     
     p <- p + geom_line(data = motif_df, aes(x = x, y = y, color = "Motif"), size = 1.5)
   }
   
   # Display the plot
   print(p)
  
  # create list with 
  return(res)

}

mot_len <- 100
mot_order <-  norder
mot_weights <- rbeta(300, 1, 3) #runif(300, -15, 15)

# Generate multiple curves with a motif structure imposed
n_curves <- 20
# general motif structure: motif 1 is in curve 1 at break 10
motif_str <- rbind.data.frame(c(1, 1, 10),
                              c(1, 1, 35),
                              c(1, 1, 89),
                              c(1, 2, 15),
                              c(1, 2, 30),
                              c(1, 4, 5),
                              c(1, 4, 60),
                              c(1, 5, 89),
                              c(1, 7, 70),
                              c(1,10, 18),
                              c(1,10, 32),
                              c(1,12, 40),
                              c(1,12, 21),
                              c(1,15, 17),
                              c(1,17, 81),
                              c(1,17, 20),
                              c(1,18, 3),
                              c(1,20, 20),
                              c(1,20, 69))
                              #c(2,17,1)
                              #c(2,17,55))
names(motif_str) <- c("motif_id", "curve", "start_break")
dim(motif_str)

# details of all the motifs to add 
mot1 <- list("mot_len" = mot_len, #length
             "dist_knots" = dist_knots, #distance between knots
             "mot_order" = mot_order, # order of the bspline
             "weight" = rbeta(300, 1, 3), # weights for the motif
             "appearance" = motif_str %>% filter(motif_id == 1)) #pattern

mot2 <- list("mot_len" = mot_len,
             "dist_knots" = dist_knots,
             "mot_order" = mot_order,
             "weight" = rbeta(300, 1, 4),
             "appearance" = motif_str %>% filter(motif_id == 2))

# list of all motifs (in this case one but it should work with more mot at the same time)
mot_details <- list(mot1) 
# mot_details <- list(mot1, mot2) 

# generate background curves
backgrounds <- lapply(1:n_curves, function(x){
  weights <- rbeta(300, 1, 3)
  generate_background_curve(len, 10, 3, weights, add_noise = TRUE)
})

# check the plot
# noise_coeff <- rnorm(length(start_point:end_point),
#                      mean = 0,
#                      sd = 0.1)
# 
# plot(backgrounds[[1]]$no_error_y, type='l')
# lines(backgrounds[[2]]$no_error_y, type='l')
# 
# ciao <- backgrounds[[1]]$no_error_y
# ciao[start_point:end_point] <- ciao[start_point:end_point] + noise_coeff
# lines(start_point:end_point,
#       ciao[start_point:end_point], type='l', col= 'red')

# --- adding the motifs
# error structure for the motifs
error_str <- rbind(
  rep(2, 110), # constant and identical
  c(rep(0.1,55), rep(2, 55)), # sd 0.1 first, sd 1 later
  c(rep(2, 55), rep(0.1, 55)), # sd 1 first, 0.1 later
  c(seq(2, 0.1, len = 55), rep(0.1, 55)) # decreasing first, constant later
)

final_curves <- backgrounds
for(i in 1:length(mot_details)){
  
  print(paste("Dealing with motif", i))
  temp_mot <- mot_details[[i]]
  curve_ids <- unique(temp_mot$appearance$curve)
  
  for(j in curve_ids){
    print(paste(" --- Adding motif", i, "to curve", j))
    temp_curve <- final_curves[[j]]
    temp_pattern <- (temp_mot$appearance %>% filter(curve == j)) %>%
      select(motif_id, start_break)
    final_curves[[j]] <- add_motif(base_curve  = temp_curve,
                                   mot_pattern = temp_pattern,
                                   mot_len     = temp_mot$mot_len,
                                   dist_knots  = temp_mot$dist_knots,
                                   mot_order   = temp_mot$mot_order,
                                   mot_weights = temp_mot$weight,
                                   error_str   = error_str
                                  )
  }
}

# compare visually the results
plot(final_curves[[1]]$no_error$no_error_y, type='l', lwd = 2)
lines(final_curves[[1]]$with_error$error_y[[1]], type = 'l', col = "grey60")
lines(final_curves[[1]]$with_error$error_y[[2]], type = 'l', col = "blue")
lines(final_curves[[1]]$with_error$error_y[[3]], type = 'l', col = "red")

# plot verification - EXPAND TO SHOW THE DIFFERENT ERROR VERSION
final_curves[[17]] #let's verify curve 17
ciao <- final_curves[[17]]$no_error$no_error_y
plot(ciao, type='l') # plot the curve
# add the motifs following the info saved in motif_info
for(q in 1:nrow(final_curves[[17]]$no_error$motif_info)){
  temp_mot <- final_curves[[17]]$no_error$motif_info[q,]
  lines(temp_mot$start_point:(temp_mot$end_point-20),
        ciao[temp_mot$start_point:(temp_mot$end_point-20)],
        type='l',
        col= temp_mot$motif_id + 1)
}

# ALL MOTIFS
my_mot_str <- lapply(final_curves,
                      function(x){
                        temp <- as.list(x)
                        print(class(temp$no_error$motif_info)) 
                        temp$no_error$motif_info
                      })

my_mot_str <- do.call(rbind, my_mot_str)

# the lines in a matrix
Y_sd0 <- lapply(final_curves,
                function(x){
                  x$no_error_y
                })
Y_sd0 <- do.call(cbind, Y_sd0)
matplot((Y_sd0), type='l')



