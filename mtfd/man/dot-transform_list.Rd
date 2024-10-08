% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4_utilities.R
\name{.transform_list}
\alias{.transform_list}
\title{Transform List Structure}
\usage{
.transform_list(lst, noise_str)
}
\arguments{
\item{lst}{A list containing elements to be transformed.}

\item{noise_str}{A string indicating the noise structure to apply.}
}
\value{
The transformed list with updated structures.
}
\description{
This function modifies a list of structures by checking for the presence of a "with_noise" field. 
If absent, it creates a new sub-list containing relevant fields.
}
\details{
Transform List Structure
}
