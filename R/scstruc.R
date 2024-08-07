#' @title scstruc
#' 
#' @description
#' The main function for scstruc package.
#' 
#' 
#' @details
#' if `Hurdle` algorithm is specified, it returns the list of `hurdle`, `undir`, `bn`.
#' `hurdle`: returned object of HurdleNormal::fitHurdle.
#' `undir`: Undirected network with the best BIC.
#' `bn`: bn object
#' 
#' @param spe SingleCellExperiment, SpatialExpetiment object
#' @param candidateGenes candidate genes to be used.
#' @param label label column name
#' @param labelName label name used to subset the cells
#' @param algorithm inference algorithm of BN
#' If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
#' If you use CCDr algorithm, the {bn} object is returned.
#' @param algorithim.args list of algorithm arguments
#' @param clusterLabel (experimental) aggregate the count across this cluster to speed up computation
#' @param useAssay use assay, default to logcounts
#' @param input default to NULL, if specified use this input matrix for inference.
#' @param rmNeg remove genes with all negative abundances
#' @param nonzero parameter to filter the genes with many zeros.
#' By default perform no filtering.
#' @param verbose used for logging.
#' @param label label column
#' @param changeSymbol automatically change the node name of input matrix to gene symbol
#' @param symbolColumn if changeSymbol is TRUE, use `symbolColumn` in rowData of the input object.
#' @param returBn By default, returns the bn class object
#' @param returnData By default, return the data used for the inference.
#' @param barcodeColumn barcode for cell identification
#' @export
scstruc <- function(spe, candidateGenes, label=NULL,
    labelName=NULL, algorithm="mmhc",
	algorithm.args=list(), returnBn=TRUE,
    returnData=TRUE, clusterLabel=NULL,
    verbose=FALSE, useAssay="logcounts",
    barcodeColumn="row", changeSymbol=TRUE,
    symbolColumn="Symbol", input=NULL,
    nonzero=1, rmNeg=FALSE) {
    ## .getInput for global label
    if (is.null(input)) {
        input <- .getInput(spe, candidateGenes, label, labelName, useAssay, barcodeColumn, clusterLabel, verbose,
            changeSymbol, symbolColumn, nonzero, rmNeg)
    }
    global_graph <- .getStruc(input, algorithm, algorithm.args, verbose)

    if (algorithm %in% c("ccdr","ccdr.boot","Hurdle")) {
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