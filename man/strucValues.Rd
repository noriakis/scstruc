% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/globalStrucValues.R
\name{strucValues}
\alias{strucValues}
\title{strucValues}
\usage{
strucValues(
  spe,
  global_tbl_graph = NULL,
  labels,
  exclude_label = NA,
  summarize_func = mean,
  variation_func = sd,
  bn = NULL,
  barcode = "barcode_id",
  assay = "logcounts",
  verbose = FALSE
)
}
\description{
Fit the network parameters using SingleCellExperiment object and bn object.
}
