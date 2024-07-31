#' scstruc
#' 
#' The main function for scstruc package.
#' 
#' @param spe SingleCellExperiment, SpatialExpetiment object
#' @param algorithm inference algorithm of BN
#' If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
#' If you use CCDr algorithm, the {bn} object is returned.
#' @param clusterLabels (experimental) aggregate the count across this cluster to speed up computation
#' @param useAssay use assay, default to logcounts
#' @param input default to NULL
#' @param label label column
#' @param returBn By default, returns the bn class object
#' @param returnData By default, return the data used for the inference.
#' @param barcodeColumn barcode for cell identification
#' @export
scstruc <- function(spe, candidateGenes, label=NULL, labelName=NULL, algorithm="mmhc",
	algorithm.args=list(), returnBn=TRUE, returnData=TRUE, clusterLabel=NULL,verbose=FALSE, useAssay="logcounts",
    barcodeColumn="row", changeSymbol=TRUE, symbolColumn="Symbol", input=NULL, nonzero=1, rmNeg=FALSE) {
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