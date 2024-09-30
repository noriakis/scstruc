#' @title scstruc
#' 
#' @description
#' The main function for scstruc package.
#' Perform input subset, structure learning.
#' 
#' @details
#' \code{algorithm}: inference algorithm (two-stage or one-stage)
#' if \code{Hurdle} algorithm is specified, it returns the list of \code{hurdle}, \code{undir}, \code{bn}.
#' \code{hurdle}: returned object of HurdleNormal::fitHurdle.
#' \code{undir}: Undirected network with the best (lowest) BIC.
#' \code{bn}: bn object
#' if \code{ccdr}, the list of network for multiple lambdas will be returned
#' 
#' @param spe SingleCellExperiment, SpatialExpetiment object
#' @param candidateGenes candidate genes to be used.
#' @param label label column name
#' @param labelName label name used to subset the cells
#' @param algorithm inference algorithm of BN
#' If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
#' If you use CCDr algorithm, the \code{bn} object is returned.
#' @param algorithm.args list of algorithm arguments
#' @param clusterLabel (experimental) aggregate the count across this cluster to speed up computation
#' @param useAssay use assay, default to logcounts
#' @param input default to NULL, if specified use this input matrix for inference.
#' @param rmNeg remove genes with all negative abundances
#' @param zeroFilt parameter to filter the genes with many zeros.
#' Genes with the sum of zero expression cells <= the number of all cells * `zeroFilt`
#' will be used as input. By default perform no filtering (1).
#' This filter will be performed *after* the genes and cells subset by the other parameters.
#' @param verbose used for logging.
#' @param label label column
#' @param changeSymbol automatically change the node name of input matrix to gene symbol
#' @param symbolColumn if changeSymbol is TRUE, use `symbolColumn` in rowData of the input object.
#' @param returnBn By default, returns the bn class object
#' @param returnData By default, return the data used for the inference.
#' @param returnOnlyData return only the data, not performing inference.
#' @importFrom dplyr %>% mutate
#' @importFrom tidygraph %E>% %N>%
#' @importFrom methods is
#' @import scran
#' @importFrom stats BIC as.formula coef coefficients
#' @importFrom stats gaussian nobs sd setNames summary.glm
#' @examples
#' library(scran)
#' library(scstruc)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' included_genes <- sample(row.names(sce), 10)
#' gs <- scstruc(sce, included_genes, changeSymbol=FALSE)
#' @export
scstruc <- function(spe, candidateGenes, label=NULL,
    labelName=NULL, algorithm="mmhc",
	algorithm.args=list(), returnBn=TRUE,
    returnData=TRUE, clusterLabel=NULL,
    verbose=FALSE, useAssay="logcounts",
    changeSymbol=TRUE, symbolColumn="Symbol", input=NULL,
    zeroFilt=1, rmNeg=FALSE, returnOnlyData=FALSE) {
    ## .getInput for global label
    if (is.null(label) & !is.null(labelName)) {
        cat_subtle("label column name not specified, defaulting to `label`\n")
        label <- "label"
    }
    if (is.null(input)) {
        input <- .getInput(spe, candidateGenes, label, labelName, useAssay, clusterLabel, verbose,
            changeSymbol, symbolColumn, zeroFilt, rmNeg)
    }
    if (returnOnlyData) {
        return(input)
    }
    global_graph <- .getStruc(input, algorithm, algorithm.args, verbose)

    if (algorithm %in% c("ccdr","ccdr.boot","Hurdle","Hurdle.boot")) {
        if (returnData) {
            return(list("net"=global_graph, "data"=input))
        } else {
            return(global_graph)
        }
    }
    if (returnBn) {
        if (returnData) {
            return(list("net"=global_graph, "data"=input))
        } else {
            return(global_graph)
        }
    }
    ## When not returning bn object, return tbl_graph
    global_tbl_graph <- global_graph |> 
        bnlearn::as.igraph() |>
        tidygraph::as_tbl_graph()
    return(global_tbl_graph)
}