% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strucValues.R
\name{strucValues}
\alias{strucValues}
\title{strucValues}
\usage{
strucValues(
  spe,
  labels,
  bn = NULL,
  net = NULL,
  exclude_label = NA,
  summarize_func = mean,
  variation_func = sd,
  barcode = "barcode_id",
  assay = "logcounts",
  verbose = FALSE,
  fit.hurdle = FALSE
)
}
\arguments{
\item{spe}{SingleCellExperiment object}

\item{labels}{label in colData. e.g., if `label` corresponding to cell cluster is specified,
the parameters will be fitted per cell basis. Multiple labels can be specified like
c("label", "SampleID") for `bn` object.}

\item{bn}{bn object}

\item{net}{tbl_graph object}

\item{exclude_label}{if specified, this label will be omitted from the fitting}

\item{summarize_func}{if `net` is provided, the summarization of subnetwork is performed by this function.}

\item{variation_func}{if `net` is provided, the summarization of subnetwork is performed by this function.}

\item{barcode}{barcode ID for cell}

\item{assay}{use assay}

\item{verbose}{control verbosity}

\item{fit.hurdle}{fit the parameter based on hurdle model}
}
\description{
Fit the network parameters using SingleCellExperiment object and bn object.
}
