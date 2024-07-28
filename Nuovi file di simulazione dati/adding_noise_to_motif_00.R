
# plot background
plot(backgrounds[[1]]$no_error_y, type='l')
#lines(backgrounds[[2]]$no_error_y, type='l')

my_background <- backgrounds[[10]]
ciao <- my_background$no_error_y
start_point <- 550
end_point <- 660
motif <- ciao[start_point:end_point]

# add noise and watch the effect
pdf("example of adding noise.pdf")
par(mfrow = c(2,2))
perc_noise <- c(0.05, 0.1, 0.2, 0.3, 0.5, 1, 1.5, 2) # noise as a percentage
for(perc_noise_i in perc_noise){
  
  print(perc_noise_i)
  
  # --- additive noise
  # noise to be added as a percentage of sd(motif) with mean 0
  noise_coeff <- rnorm(length(start_point:end_point),
                     mean = 0,
                     sd = sd(motif)*perc_noise_i)
  # add noise
  prova_add_err <- motif + noise_coeff
  
  # verify with a plot
  #plot(motif, type='l')
  #lines(prova_add_err, col='red')
  
  # smooth it again 
  mot_breaks <- seq(from = 1, length(prova_add_err), by = dist_knots)
  basisobj = create.bspline.basis(norder = 3, breaks = mot_breaks)
  #plot(basisobj)
  ys = smooth.basis(argvals = 1:(length(prova_add_err)),
                  y = prova_add_err,
                  fdParobj = basisobj)
  # smoothed again after adding error
  yy <- eval.fd(1:(length(prova_add_err)), ys$fd)
  #lines(yy, col='blue')
  
  # --- weighted noise
  #motif_info <- backgrounds[[1]]$motif_info
  curve_weights <- my_background$or_coeff # bspline weight (in general)
  mot_weights <- curve_weights[start_break:end_break] # bspline weight for the motif

  #  ---  error added to the weights
  noise_coeff <- rnorm(length(mot_weights),
                       mean = 0,
                       sd = rep((mot_weights)*perc_noise_i, length(coeff)))#sd(mot_weights)*perc_noise_i)
  coeff <- mot_weights + noise_coeff
  curve_weights[start_break:end_break] <- coeff
  # create fd object
  fd_curve <- fd(curve_weights, my_background$basis)
  # create the curve vector with no noise
  prova_w_err <- generate_curve_vector(fd_curve = fd_curve)
  prova_w_err <- prova_w_err[start_point:end_point]
  
  # original plot
  plot(my_background$no_error_y[start_point:end_point], type='l', lwd = 3,
       main = paste("SNR: ", var(motif)/var(yy - my_background$no_error_y[start_point:end_point]))) #main = paste("Level Noise: ", perc_noise_i*100)) # original
  # additive error
  lines(yy, type='l', col = 'red') # smooth
  lines(prova_add_err, type='l', col="grey60")
  # coeff error
  #lines(prova_w_err, type='l', col = "blue", main = paste("Level Noise: ", perc_noise_i*100, "% of Motif Sdev"))
}

dev.off()

for(i in 1:4){
  plot(prova_add_err, type='l', main = paste("SNR:", perc_noise_i))
}




# How to measure?
# - signal-to-noise ratio

# stima dell'errore post smoothing
exp_epsilon <- yy - motif
# Calcolare il SNR: varianza segnale vero (motivo pulito) diviso
# varianza stima dell'errore (motivo sporcato - motivo pulito)
SNR <- var(motif) / var(exp_epsilon)

# Convertire il SNR in decibel
SNR_dB <- 10 * log10(SNR)

SNR
SNR_dB


