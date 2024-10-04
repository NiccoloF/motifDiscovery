library(mtfd)
########################## RUN PROBKMA ######################################### 

diss = 'd0_d1_L2' 
alpha = 0.5
max_gap = 0 # no gaps allowed
iter4elong = 100 # perform elongation
trials_elong = 201 # try all possible elongations
c_max = 71 # maximum motif length 70
#run probKMA multiple times (2x3x10=60 times)
K = c(2,3) # number of clusters to try #c(2,3)
c = c(61,51) # minimum motif lengths to try #c(61,51)
n_init = 10 # number of random initializations to try

data("simulated200")

find_candidate_motifs_results = mtfd::clusterMotif(Y0=simulated200$Y0,method="ProbKMA",stopCriterion="max",
                                                   name = './results_ProbKMA_VectorData/',plot = TRUE,
                                                   probKMA_options = list(Y1=simulated200$Y1,K=K, c=c,n_init=n_init,
                                                                           names_var = 'x(t)',diss=diss,alpha=alpha),
                                                   worker_number = NULL)

################################################################################




########################## RUN FUNBIALIGN ######################################

funBialignResult = mtfd::clusterMotif(Y0=simulated200$Y0,method="FunBIalign",stopCriterion = 'fMRS',
                                      name = './results_FunBialign',plot = TRUE,
                                      funBIalign_options = list(portion_len=60,min_card=3,cut_off = NULL))
funBialignResult = mtfd::clusterMotif(Y0=simulated200$Y0,method="FunBIalign",stopCriterion = 'variance',
                                      name = './results_FunBialign',plot = TRUE,
                                      funBIalign_options = list(portion_len=60,min_card=3,cut_off = 10))

################################################################################



########################## RUN PROBKMA MATRIX ######################################

diss='d0_L2'
alpha=0
max_gap = 0.2 # no gaps allowed
iter4elong = 100 # perform elongation
trials_elong = 100
c_max = 150 
### run probKMA multiple times (2x3x10=60 times)
K = c(2, 3) # number of clusters to try
c = c(40, 50) # minimum motif lengths to try
n_init = 10 # number of random initializations to try

load("./Y.RData")

find_candidate_motifs_results = mtfd::clusterMotif(Y0=Y$Y0,method="ProbKMA",stopCriterion="max",
                                                   name = './results_ProbKMA_MatrixData/',plot = TRUE,
                                                   probKMA_options = list(K=K,c=c,n_init=n_init,diss=diss,alpha=alpha,
                                                                          iter_max = 250,c_max = c_max,iter4elong = iter4elong,
                                                                          trials_elong = trials_elong,max_gap = max_gap),
                                                   worker_number = NULL)

################################################################################
############# 
### CASE 0: NO MOTIFS
mot_len <- 100
mot_details <- NULL # or list()

builder <- mtfd::motifSimulationBuilder(N = 20,len = 300,mot_details)

curves <- mtfd::generateCurves(builder) 

mtfd::plot_motifs(builder,curves,"plots")

result <- mtfd::to_motifDiscovery(curves)
rm(list = ls())



### CASE 1: SET THE POSITION AND POINTWISE ERROR
mot_len <- 100
motif_str <- rbind.data.frame(c(1, 1, 20),
                              c(2, 1, 2), #9 
                              c(1, 3, 1),
                              c(1, 2, 1),
                              c(1, 2, 15),
                              c(1, 4, 1),
                              c(2, 5, 1),
                              c(2, 7, 1),
                              c(2,17,1))

names(motif_str) <- c("motif_id", "curve","start_break_pos")
dim(motif_str)
# cambiare da pesi a coefficienti 
# cambiare posizione nodi a posizione assoluta
# appearance --> occurrences
mot1 <- list("len" = mot_len, #length
             "coeffs" = NULL, # weights for the motif
             "occurrences" = motif_str %>% filter(motif_id == 1)) #pattern motif_str %>% filter(motif_id == 1)

mot2 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = motif_str %>% filter(motif_id == 2))

mot_details <- list(mot1,mot2)

# VECTOR ERROR 
noise_str <- c(0.1,0.5,1.0)
noise_str <- matrix(c(0.1,0.5,1.0,5.0),ncol = 1)

# MATRIX ERROR 
noise_str <- list(rbind(
  rep(2, 100), # constant and identical
  c(rep(0.1,50), rep(2, 50)), # sd 0.1 first, sd 1 later
  c(rep(2, 50), rep(0.1, 50)), # sd 1 first, 0.1 later
  c(seq(2, 0.1, len = 50), rep(0.1, 50))),
  rbind(
    rep(0.0, 100), 
    rep(0.5, 100), 
    rep(1.0, 100), 
    rep(5.0, 100)))
  
builder <- mtfd::motifSimulationBuilder(N = 20,len = 300,mot_details,
                                        distribution = 'beta')

curves <- mtfd::generateCurves(builder,noise_type = 'pointwise', noise_str = noise_str,seed_motif = 321) 

mtfd::plot_motifs(builder,curves,name="plots")

result <- mtfd::to_motifDiscovery(curves)

rm(list = ls())

### CASE 2: SET THE POSITION AND COEFF ERROR
mot_len <- 100
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
# cambiare da pesi a coefficienti 
# cambiare posizione nodi a posizione assoluta
# appearance --> occurrences
mot1 <- list("len" = mot_len, #length
             "coeffs" = NULL, # weights for the motif
             "occurrences" = motif_str %>% filter(motif_id == 1)) #pattern motif_str %>% filter(motif_id == 1)

mot2 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = motif_str %>% filter(motif_id == 2))

mot_details <- list(mot1,mot2)

# SCALAR ERROR
noise_str <- 0.1

# VECTOR ERROR 
noise_str <- c(0.1,1.0,5.0)
noise_str <- list(c(0.1,5.0,10.0),c(0.0,0.0,0.0))
noise_str <- matrix(c(0.1,0.5,1.0,5.0),ncol = 1)

builder <- mtfd::motifSimulationBuilder(N = 20,len = 300,mot_details,
                                        distribution = 'beta')

curves <- mtfd::generateCurves(builder,noise_type = 'coeff',noise_str,only_der = FALSE)

mtfd::plot_motifs(builder,curves,"plots")

result <- mtfd::to_motifDiscovery(curves)

### CASE 3: RANDOM POSITION AND POINTWISE ERROR
mot_len <- 100

mot1 <- list("len" = mot_len, #length
             "coeffs" = NULL, # weights for the motif
             "occurrences" = 5) #pattern motif_str %>% filter(motif_id == 1)

mot2 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = 6)
mot3 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = 4)

mot_details <- list(mot1,mot2)

# SCALAR ERROR
noise_str <- 0.5

# VECTOR ERROR 
noise_str <- c(0.1,0.5,1.0)
noise_str <- matrix(c(0.1,0.5,1.0,5.0),ncol = 1)

builder <- mtfd::motifSimulationBuilder(N = 100,len = 300,mot_details,
                                        distribution = 'beta')

curves <- mtfd::generateCurves(builder,noise_type = 'pointwise',noise_str)

mtfd::plot(builder,curves,getwd())



### CASE 4: RANDOM POSITION AND COEFF ERROR
mot_len <- 100
mot1 <- list("len" = mot_len, #length
             "weights" = NULL, # weights for the motif
             "occurrences" = 5) #pattern motif_str %>% filter(motif_id == 1)

mot2 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = 6)
mot3 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = 4)

mot_details <- list(mot1,mot2)

# SCALAR ERROR
noise_str <- 0.1

# VECTOR ERROR 
noise_str <- list(c(0.0,1.0,5.0),c(0.1,2.0,4.0))
noise_str <- matrix(c(0.1,0.5,1.0,5.0),ncol = 1)

builder <- mtfd::motifSimulationBuilder(N = 20,len = 300,mot_details,
                                        distribution = 'beta')

curves <- mtfd::generateCurves(builder,noise_type = 'coeff',noise_str)

mtfd::plot_motifs(builder,curves,"plots")















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
# cambiare da pesi a coefficienti 
# cambiare posizione nodi a posizione assoluta
# appearance --> occurrences
mot1 <- list("len" = mot_len, #length
             "weights" = NULL, # weights for the motif
             "appearance" = motif_str %>% filter(motif_id == 1)) #pattern motif_str %>% filter(motif_id == 1)

mot2 <- list("len" = mot_len,
             "weights" = NULL,
             "appearance" = motif_str %>% filter(motif_id == 2))
mot3 <- list("len" = mot_len,
             "weights" = NULL,
             "appearance" = 4)

mot_details <- list(mot1,mot2)
#############   

#############   
noise_str <- list(rbind(rep(2, 100),c(rep(0.1,50),rep(2, 50))), 
                  rbind(rep(0.0, 100),rep(5.0, 100)))

##############  


# Mettere min_dist_motifs a NULL 
# norder = 3 
builder <- mtfd::motifSimulationBuilder(N = 20,len = 300,mot_details,
                                        distribution = 'beta')
# aggiungere flag sulla tipologia di errore da aggiungere 
# cambiare no_error -> background
# aggiungere curva con motif ma no errore -> no_noise
# with_noise -> curva con motivi + errore
# aggiungere SNR su singolo motivo 
curves <- mtfd::generateCurves(builder,noise_type = 'pointwise',noise_str)
curves <- mtfd::generateCurves(builder,noise_type = 'coeff',noise_str)
# Aggiungere errori singoli per motivi 

mtfd::plot_motifs(builder,curves,"plots")

mtfd::motifSimulationApp(noise_str,mot_details)

# Aggiungere plots occorrenze dei motivi con il motivo senza errore.

# PRIMA DI FAR GIRARE FIND CANDIDATE SPEZZO LA CURVA IN MODO RANDOM LA CURVA IN
# IL MOTIF SEARCH VA FATTO SU TUTTA LA CURVA E NON QUELLA SPEZZETTATO 
# METTERE ALL'INTERNO DEI CHECK PER CAPIRE SE IL NUMERO DI CURVE DA GENERARE è SUFFICIENTI
mot_len <- 100
motif_str <- rbind.data.frame(c(1, 1,80),
                              c(1,1,1),
                              c(2,1,30),
                              c(2,1,60))

names(motif_str) <- c("motif_id", "curve","start_break_pos")
dim(motif_str)
# cambiare da pesi a coefficienti 
# cambiare posizione nodi a posizione assoluta
# appearance --> occurrences
mot1 <- list("len" = mot_len, #length
             "coeffs" = NULL, # weights for the motif
             "occurrences" = motif_str %>% filter(motif_id == 1)) #pattern motif_str %>% filter(motif_id == 1)
mot2 <- list("len" = mot_len, #length
             "coeffs" = NULL, # weights for the motif
             "occurrences" = motif_str %>% filter(motif_id == 2))

mot_details <- list(mot1,mot2)

noise_str <- list(c(3.0),c(1.0))

builder <- mtfd::motifSimulationBuilder(N=1,len = 1000,mot_details = mot_details,norder = 3,
                                        coeff_min=-15,coeff_max=15,
                                        dist_knots=10,min_dist_motifs=NULL,
                                        distribution = 'beta')

curves <- mtfd::generateCurves(builder,noise_type = 'coeff', noise_str = noise_str,seed_motif = 100,only_der = FALSE) 

mtfd::plot_motifs(builder,curves,"test")

Y0 <- mtfd::to_motifDiscovery(curves)

diss='d0_L2'
alpha=0
max_gap = 0.2 
iter4elong = 100 # perform elongation
trials_elong = 100
c_max = 120
### run probKMA multiple times (2x3x10=60 times)
K = c(2,3) # number of clusters to try
c = c(80,100) # minimum motif lengths to try
n_init = 10 # number of random initializations to try

find_candidate_motifs_results = mtfd::clusterMotif(Y0=Y0[[1]],method="ProbKMA",stopCriterion="max",
                                                   name = './results_ProbKMA_singleData/',plot = TRUE,
                                                   probKMA_options = list(ter4clean=50,K=K, c=c,n_init=n_init,
                                                                          names_var = 'x(t)',diss=diss,alpha=alpha,
                                                                          sil_threshold=0.8,n_subcurves=5),
                                                   worker_number = NULL)




