% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scstruc.R
\name{scstruc}
\alias{scstruc}
\title{scstruc}
\usage{
scstruc(
  spe,
  candidateGenes,
  label = NULL,
  labelName = NULL,
  algorithm = "mmhc",
  algorithm.args = list(),
  returnBn = TRUE,
  returnData = TRUE,
  clusterLabel = NULL,
  verbose = FALSE,
  useAssay = "logcounts",
  barcodeColumn = "row",
  changeSymbol = TRUE,
  symbolColumn = "Symbol",
  input = NULL,
  nonzero = 1,
  rmNeg = FALSE
)
}
\arguments{
\item{spe}{SingleCellExperiment, SpatialExpetiment object}

\item{candidateGenes}{candidate genes to be used.}

\item{label}{label column}

\item{labelName}{label name used to subset the cells}

\item{algorithm}{inference algorithm of BN
If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
If you use CCDr algorithm, the {bn} object is returned.}

\item{returnData}{By default, return the data used for the inference.}

\item{clusterLabel}{(experimental) aggregate the count across this cluster to speed up computation}

\item{verbose}{used for logging.}

\item{useAssay}{use assay, default to logcounts}

\item{barcodeColumn}{barcode for cell identification}

\item{changeSymbol}{automatically change the node name of input matrix to gene symbol}

\item{symbolColumn}{if changeSymbol is TRUE, use `symbolColumn` in rowData of the input object.}

\item{input}{default to NULL, if specified use this input matrix for inference.}

\item{nonzero}{parameter to filter the genes with many zeros.
By default perform no filtering.}

\item{rmNeg}{remove genes with all negative abundances}

\item{algorithim.args}{list of algorithm arguments}

\item{returBn}{By default, returns the bn class object}
}
\description{
The main function for scstruc package.
}
\details{
if `Hurdle` algorithm is specified, it returns the list of `hurdle`, `undir`, `bn`.
`hurdle`: returned object of HurdleNormal::fitHurdle.
`undir`: Undirected network with the best BIC.
`bn`: bn object
}
