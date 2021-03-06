% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dyadmlvar.R
\name{lasso}
\alias{lasso}
\title{Use a LASSO algorithm to find important predictors}
\usage{
lasso(Data, Predictors, Outcome, Seeds = 1, Train = F,
  PropOfTrain = 0.75, Plot = F)
}
\arguments{
\item{Data}{a dataframe with columns for each predictor and for the outcome
variable}

\item{Predictors}{a list of strings with predictor variables names}

\item{Outcome}{a string with the name of an outcome variable}

\item{Seeds}{a seed for randomly selecting some part of the data to be
training data and some to be testing data. default -}

\item{Train}{Should we split data to training and test datasets (FALSE uses
all data for both)}

\item{PropOfTrain}{How much of the data to use for training the model}

\item{Plot}{shoud we show a plot of the parameter number and log likelihood?}
}
\value{
a list of estimates for the effect of the valid predictors, and a R^2
  statistic for the final model
}
\description{
\code{lasso} Function that uses the LASSO shrinkage reduction method
to find meaningful predictors of an outcome vector
}
\examples{
ntwrk <- read_network(fit1)
#intra-partner variables from the partner B
intra_B_vars <- get_names(ntwrk, part = INTRA_B, time = ALL_NETWORK)
full_data <- merge(ntwrk$all, sat, by = "ID", all=TRUE, sort=TRUE)
lasso(full_data, intra_B_vars, "W_csi_resid")
}
