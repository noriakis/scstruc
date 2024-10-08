% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccdr.boot.R
\name{ccdr.boot}
\alias{ccdr.boot}
\title{ccdr.boot}
\usage{
ccdr.boot(
  data,
  R = 200,
  m = nrow(data),
  lambdas.length = 20,
  alpha = 2,
  gamma = 2,
  return.all = FALSE
)
}
\arguments{
\item{data}{data (row: sample, column: gene)}

\item{R}{replicate number}

\item{m}{sampling number}

\item{lambdas.length}{passed to generate.lambdas}

\item{alpha}{passed to ccdr.run, default to 2}

\item{gamma}{passed to ccdr.run, default to 2}

\item{return.all}{return all the network}
}
\value{
list of the bn.strength per lambda
}
\description{
Bootstrap-based arc strength calculation based on CCDr algorithm.
The function uses ccdrAlgorithm::ccdr.run for the sampled data, and returns the 
bootstrap strength for the arcs per lambda. The same lambdas are used for all the inference.
The networks with the same lambda are used to compute bootstapped network.
}
