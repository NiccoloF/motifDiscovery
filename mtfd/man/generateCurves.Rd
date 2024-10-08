% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_methods.R
\name{generateCurves}
\alias{generateCurves}
\title{Generate Functional Curves with Embedded Motifs}
\usage{
generateCurves(
  object,
  noise_type = NULL,
  noise_str = NULL,
  seed_background = 777,
  seed_motif = 43213,
  only_der = TRUE,
  coeff_min_shift = -10,
  coeff_max_shift = 10
)
}
\arguments{
\item{object}{An S4 object of class `motifSimulation` that has been previously constructed using the `motifSimulationBuilder` function. This object encapsulates all necessary parameters and configurations for curve and motif generation (mandatory).}

\item{noise_type}{A character string specifying the type of noise to add to the curves. Acceptable values are `'pointwise'` for adding noise directly to data points or `'coeff'` for perturbing the coefficients of the basis functions (mandatory).}

\item{noise_str}{A list detailing the structure and magnitude of the noise to be added for each motif. 
\itemize{
  \item If `noise_type` is `'pointwise'`, `noise_str` should contain vectors or matrices indicating the noise level for each motif.
  \item If `noise_type` is `'coeff'`, `noise_str` should include individual values or vectors representing the noise to be applied to the coefficients.
}
This parameter allows fine-grained control over the noise characteristics applied to each motif (mandatory).}

\item{seed_background}{An integer value setting the seed for the random number generator used in background curve generation. This ensures reproducibility of the background curves. Default is `777`.}

\item{seed_motif}{An integer value setting the seed for the random number generator used in motif generation. This ensures reproducibility of the motif embedding process. Default is `43213`.}

\item{only_der}{A logical value indicating whether to apply only derivative-based modifications to the motifs (`TRUE`) or to add a vertical shift in addition to derivative modifications (`FALSE`). Setting to `FALSE` introduces vertical shifts to each motif instance, allowing motifs to appear at different baseline levels within the curves. Default is `TRUE`.}

\item{coeff_min_shift}{A numeric value specifying the minimum vertical shift to be applied to motifs when `only_der` is set to `FALSE`. This parameter controls the lower bound of the vertical displacement of motifs. Default is `-10`.}

\item{coeff_max_shift}{A numeric value specifying the maximum vertical shift to be applied to motifs when `only_der` is set to `FALSE`. This parameter controls the upper bound of the vertical displacement of motifs. Default is `10`.}
}
\value{
A list containing the following components:
\itemize{
  \item \strong{basis}: The basis functions used to represent the curves.
  \item \strong{background}: A list containing the coefficients and the background curves without any motifs.
  \item \strong{no_noise}: A list containing the coefficients and the background curves with embedded motifs but without added noise.
  \item \strong{with_noise}: A list containing the noise structure and the curves with embedded motifs and added noise.
  \item \strong{SNR}: A list of Signal-to-Noise Ratio (SNR) metrics calculated for each motif within each curve, useful for assessing the quality of motif embedding.
}
}
\description{
The `generateCurves` function is designed to create synthetic functional data by embedding predefined motifs into background curves. This is particularly useful for testing and benchmarking motif discovery and clustering algorithms in functional data analysis. By allowing the incorporation of various noise types and controlled motif placements, the function provides a flexible framework for simulating realistic scenarios where motifs may or may not be noisy.

The function supports two types of noise addition:
\itemize{
  \item \strong{Pointwise Noise}: Adds noise directly to the data points of the curves, simulating random fluctuations or measurement errors.
  \item \strong{Coefficient Noise}: Perturbs the coefficients of the basis functions used to represent the curves, allowing for smoother variations and controlled distortions.
}

Additionally, `generateCurves` allows for the specification of vertical shifts to motifs, enabling the simulation of motifs appearing at different baseline levels within the curves. The function ensures that all generated subcurves meet the minimum motif length requirement, maintaining the integrity of the embedded motifs.

This function is integral to the `mtfd` package's motif simulation capabilities, providing users with the ability to create complex functional datasets tailored to their specific research or testing needs.
}
\examples{
\dontrun{
# Example 0: Special case with no motifs
mot_len <- 100
mot_details <- NULL  # or list()
builder <- mtfd::motifSimulationBuilder(N = 20, len = 300, mot_details)
curves <- mtfd::generateCurves(builder)

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
  "occurrences" = motif_str \%>\% filter(motif_id == 1)
)
mot2 <- list(
  "len" = mot_len, 
  "coeffs" = NULL, 
  "occurrences" = motif_str \%>\% filter(motif_id == 2)
)
mot_details <- list(mot1, mot2)

# Define noise structure for pointwise noise
noise_str <- list(
  rbind(rep(2, 100), rep(c(rep(0.1, 50), rep(2, 50)), 1)),
  rbind(rep(0.0, 100), rep(0.5, 100))
)

# Build the simulation object
builder <- mtfd::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')

# Generate curves with pointwise noise
curves <- mtfd::generateCurves(builder, noise_type = 'pointwise', noise_str = noise_str)

# Example 2: Set the motif position and add coefficient noise
# Define noise structure for coefficient noise
noise_str <- list(c(0.1, 1.0, 5.0), c(0.0, 0.0, 0.0))

# Generate curves with coefficient noise without vertical shifts
curves <- mtfd::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)

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
builder <- mtfd::motifSimulationBuilder(N = 20, len = 300, mot_details, distribution = 'beta')

# Generate curves with pointwise noise and vertical shifts
curves <- mtfd::generateCurves(builder, noise_type = 'pointwise', noise_str, only_der = FALSE)

# Example 4: Random motif positions and add coefficient noise
# Define noise structure for coefficient noise
noise_str <- list(c(0.1, 5.0, 10.0), c(0.1, 5.0, 10.0))

# Generate curves with coefficient noise and vertical shifts
curves <- mtfd::generateCurves(builder, noise_type = 'coeff', noise_str, only_der = FALSE)
}
}
