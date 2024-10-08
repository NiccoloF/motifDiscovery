% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_utilities.R
\name{add_error_to_motif}
\alias{add_error_to_motif}
\title{Add Additive Error to Motif}
\usage{
add_error_to_motif(or_y, noise_str, start_point, end_point, k)
}
\arguments{
\item{or_y}{A numeric vector representing the original motif values.}

\item{noise_str}{A string or numeric vector indicating the structure of the noise to be added.}

\item{start_point}{A numeric vector indicating the starting points of motifs.}

\item{end_point}{A numeric vector indicating the ending points of motifs.}

\item{k}{An integer representing the current iteration or motif index.}
}
\value{
A numeric vector of the motif values with added noise.
}
\description{
This function adds an additive noise to a given motif based on specified parameters. The noise is generated as a percentage of the standard deviation of the motif.
}
\details{
Add Error to Motif
}
