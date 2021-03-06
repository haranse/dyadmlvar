% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dyadmlvar.R
\name{get_names}
\alias{get_names}
\title{Get a list of variable names from a dyadNetwork object}
\usage{
get_names(ntwrk, part = ALL_NETWORK, time = ALL_NETWORK)
}
\arguments{
\item{ntwrk}{a dyadNetwork object}

\item{part}{which kind of variables to retrieve:
INTRA_A - interconnections between partner A variables
INTRA_B - interconnections between partner B variables
INTER - connections between partner A to partner B variables
STRENGTH - node strengths
DENSITY - density of intra-partner A network, intra-partner B
   network and inter-partner networks and the ratio between
   inter-partner density to the average of intra-partner densities
ALL_NETWORK - all variables}

\item{time}{which network to retreive variables from:
TEMPORAL - temporal network
CONTEMP - contemporaneous network
ALL_NETWORK - both}
}
\value{
a list of variable name strings
}
\description{
\code{get_names} Get a list of variable names from a dyadNetwork object, e.g.
variable names for intra-partner connections for partner A
}
\examples{
ntwrk <- read_network(fit1)
#all variable names
all_vars <- get_names(ntwrk)
#inter-partner variables from the contemporaneous network
inter_contemp_vars <- get_names(ntwrk, part = INTER, time = CONTEMP)
}
