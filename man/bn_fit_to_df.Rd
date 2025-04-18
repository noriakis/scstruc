% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bn_fit_to_df}
\alias{bn_fit_to_df}
\title{bn_fit_to_df}
\usage{
bn_fit_to_df(fitted)
}
\arguments{
\item{fitted}{bn.fit object}
}
\description{
This function coverts the bn.fit object to data.frame.
Only the fitted parameters are retained and nodes not in the fitted object will not be included.
}
\examples{
data(gaussian.test)
test <- head(gaussian.test)
struc <- hc(test)
bn_fit_to_df(bn.fit(struc, test))
}
