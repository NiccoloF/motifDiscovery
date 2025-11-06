pkgname <- "funMoDisco"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "funMoDisco-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('funMoDisco')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("discoverMotifs")
### * discoverMotifs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: discoverMotifs
### Title: Functional Motif Discovery
### Aliases: discoverMotifs

### ** Examples

## No test: 
# Example 1: Discover motifs using ProbKMA

# Define dissimilarity measure and weight parameter
diss <- 'd0_d1_L2'
alpha <- 0.5

# Define number of motifs and their minimum lengths
K <- c(2, 3)
c <- c(61, 51)
n_init <- 3

# Load simulated data
data("simulated200")

# Perform motif discovery using ProbKMA
results <- funMoDisco::discoverMotifs(
  Y0 = simulated200$Y0,
  method = "ProbKMA",
  stopCriterion = "max",
  name = tempdir(),
  plot = FALSE,
  probKMA_options = list(
    Y1 = simulated200$Y1,
    K = K,
    c = c,
    n_init = n_init,
    diss = diss,
    alpha = alpha
  ),
  worker_number = 1
)

# Example 2: Discover motifs using funBIalign
results_funbialign <- funMoDisco::discoverMotifs(
  Y0 = simulated200$Y0,
  method = "funBIalign",
  stopCriterion = "fMSR",
  name = getwd(),
  plot = FALSE,
  funBIalign_options = list(
    portion_len = 60,
    min_card = 3,
    cut_off = 1.0
  ),
  worker_number = 1
  ) 

## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("discoverMotifs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generateCurves-motifSimulation-method")
### * generateCurves-motifSimulation-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generateCurves,motifSimulation-method
### Title: Generate Functional Curves with Embedded Motifs
### Aliases: generateCurves,motifSimulation-method

### ** Examples

## No test: 
# Example 0: Special case with no motifs
mot_len <- 100
mot_details <- NULL  # or list()
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details)
curves <- funMoDisco::generateCurves(builder)

# Example 1: Set the motif position and add pointwise noise
# Define motif positions and their respective curves
motif_str <- rbind.data.frame(
  c(1, 1, 20),
  c(2, 1, 2),
  c(2, 7, 1),
  c(2,17,1)
)
names(motif_str) <- c("motif_id", "curve", "start_break_pos")

# Define motif details
mot1 <- list(
  "len" = mot_len, 
  "coeffs" = NULL, 
  "occurrences" = motif_str %>% filter(motif_id == 1)
)
mot2 <- list(
  "len" = mot_len, 
  "coeffs" = NULL, 
  "occurrences" = motif_str %>% filter(motif_id == 2)
)
mot_details <- list(mot1, mot2)

# Define noise structure for pointwise noise
noise_str <- list(
  rbind(rep(2, 100), rep(c(rep(0.1, 50), rep(2, 50)), 1)),
  rbind(rep(0.0, 100), rep(0.5, 100))
)

# Build the simulation object
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')

# Generate curves with pointwise noise
curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str = noise_str)

# Example 2: Set the motif position and add coefficient noise
# Define noise structure for coefficient noise
noise_str <- list(c(0.1, 1.0, 5.0), c(0.0, 0.0, 0.0))

# Generate curves with coefficient noise without vertical shifts
curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)

# Example 3: Random motif positions and add pointwise noise
mot1 <- list(
  "len" = mot_len,
  "coeffs" = NULL,
  "occurrences" = 5
)
mot2 <- list(
  "len" = mot_len,
  "coeffs" = NULL,
  "occurrences" = 6
)
mot_details <- list(mot1, mot2)

# Define noise structure for pointwise noise
noise_str <- list(
  rbind(rep(2, 100)),
  rbind(rep(0.5, 100))
)

# Build the simulation object
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')

# Generate curves with pointwise noise and vertical shifts
curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str, only_der = FALSE)

# Example 4: Random motif positions and add coefficient noise
# Define noise structure for coefficient noise
noise_str <- list(c(0.1, 5.0, 10.0), c(0.1, 5.0, 10.0))

# Generate curves with coefficient noise and vertical shifts
curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generateCurves-motifSimulation-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generateCurves")
### * generateCurves

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generateCurves
### Title: Generate Functional Curves with Embedded Motifs
### Aliases: generateCurves

### ** Examples

## No test: 
# Example 0: Special case with no motifs
mot_len <- 100
mot_details <- NULL  # or list()
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details)
curves <- funMoDisco::generateCurves(builder)

# Example 1: Set the motif position and add pointwise noise
# Define motif positions and their respective curves
motif_str <- rbind.data.frame(
  c(1, 1, 20),
  c(2, 1, 2),
  c(2, 7, 1),
  c(2,17,1)
)
names(motif_str) <- c("motif_id", "curve", "start_break_pos")

# Define motif details
mot1 <- list(
  "len" = mot_len, 
  "coeffs" = NULL, 
  "occurrences" = motif_str %>% filter(motif_id == 1)
)
mot2 <- list(
  "len" = mot_len, 
  "coeffs" = NULL, 
  "occurrences" = motif_str %>% filter(motif_id == 2)
)
mot_details <- list(mot1, mot2)

# Define noise structure for pointwise noise
noise_str <- list(
  rbind(rep(2, 100), rep(c(rep(0.1, 50), rep(2, 50)), 1)),
  rbind(rep(0.0, 100), rep(0.5, 100))
)

# Build the simulation object
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')

# Generate curves with pointwise noise
curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str = noise_str)

# Example 2: Set the motif position and add coefficient noise
# Define noise structure for coefficient noise
noise_str <- list(c(0.1, 1.0, 5.0), c(0.0, 0.0, 0.0))

# Generate curves with coefficient noise without vertical shifts
curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)

# Example 3: Random motif positions and add pointwise noise
mot1 <- list(
  "len" = mot_len,
  "coeffs" = NULL,
  "occurrences" = 5
)
mot2 <- list(
  "len" = mot_len,
  "coeffs" = NULL,
  "occurrences" = 6
)
mot_details <- list(mot1, mot2)

# Define noise structure for pointwise noise
noise_str <- list(
  rbind(rep(2, 100)),
  rbind(rep(0.5, 100))
)

# Build the simulation object
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')

# Generate curves with pointwise noise and vertical shifts
curves <- funMoDisco::generateCurves(builder, noise_type = 'pointwise', noise_str, only_der = FALSE)

# Example 4: Random motif positions and add coefficient noise
# Define noise structure for coefficient noise
noise_str <- list(c(0.1, 5.0, 10.0), c(0.1, 5.0, 10.0))

# Generate curves with coefficient noise and vertical shifts
curves <- funMoDisco::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generateCurves", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("motifSimulationBuilder")
### * motifSimulationBuilder

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: motifSimulationBuilder
### Title: Create motifSimulation Object
### Aliases: motifSimulationBuilder

### ** Examples

## No test: 
mot_len <- 100
motif_str <- rbind.data.frame(c(1, 1, 20),
                              c(2, 1, 2), 
                              c(1, 3, 1),
                              c(1, 2, 1),
                              c(1, 2, 15),
                              c(1, 4, 1),
                              c(2, 5, 1),
                              c(2, 7, 1),
                              c(2,17,1))

names(motif_str) <- c("motif_id", "curve","start_break_pos")

mot1 <- list("len" = mot_len, #length
             "coeffs" = NULL, # weights for the motif
             "occurrences" = motif_str %>% filter(motif_id == 1))

mot2 <- list("len" = mot_len,
             "coeffs" = NULL,
             "occurrences" = motif_str %>% filter(motif_id == 2))

mot_details <- list(mot1,mot2)

# MATRIX ERROR 
noise_str <- list(rbind(rep(2, 100)),
                  rbind(rep(0.0, 100)))

builder <- funMoDisco::motifSimulationBuilder(N = 20,len = 300, mot_details,
                                        distribution = 'beta')
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("motifSimulationBuilder", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_motifs-motifSimulation-method")
### * plot_motifs-motifSimulation-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_motifs,motifSimulation-method
### Title: Plot Embedded Motifs in Functional Curves
### Aliases: plot_motifs,motifSimulation-method

### ** Examples

## No test: 
# Example: Plotting motifs in generated curves
# Assume 'builder' has been created and 'curves' have been generated using generateCurves
mot_details <- NULL
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details)
curves <- funMoDisco::generateCurves(builder)

# Specify the directory to save plots
plots_name <- "plots_1"  

# Generate and save the plots
funMoDisco::plot_motifs(builder, curves, plots_name, path = tempdir())
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_motifs-motifSimulation-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_motifs")
### * plot_motifs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_motifs
### Title: Plot Embedded Motifs in Functional Curves
### Aliases: plot_motifs

### ** Examples

## No test: 
# Example: Plotting motifs in generated curves
# Assume 'builder' has been created and 'curves' have been generated using generateCurves

mot_details <- NULL
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details = mot_details)
curves <- funMoDisco::generateCurves(builder)

# Specify the directory to save plots
plots_name <- "plots_1"  

# Generate and save the plots
funMoDisco::plot_motifs(builder, curves, plots_name, path = tempdir())
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_motifs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("to_motifDiscovery-list-method")
### * to_motifDiscovery-list-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: to_motifDiscovery,list-method
### Title: to_motifDiscovery
### Aliases: to_motifDiscovery,list-method

### ** Examples

## No test: 
builder <- funMoDisco::motifSimulationBuilder(N = 20, len = 300, mot_details = NULL)
curves <- funMoDisco::generateCurves(builder)
formatted_curves <- to_motifDiscovery(curves)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("to_motifDiscovery-list-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("to_motifDiscovery")
### * to_motifDiscovery

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: to_motifDiscovery
### Title: to_motifDiscovery
### Aliases: to_motifDiscovery

### ** Examples

## No test: 
N <- 20 # number of curves
len <- 300 # length of the curves
mot_details <- NULL # or list()

builder <- motifSimulationBuilder(N = N, len = len, mot_details = NULL)
curves <- funMoDisco::generateCurves(builder)
formatted_curves <- to_motifDiscovery(curves)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("to_motifDiscovery", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
