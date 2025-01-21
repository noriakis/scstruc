#' @title superCellMat
#' 
#' @description Return the metacell abundance matrix
#' Should provide SingleCellExperiment with logcounts data filled.
#' 
#' @param sce sce object
#' @param rank rank to be used in prcomp
#' @param pca use PCA or not
#' @param prop if genes are not specified, `getTopHVGs` will run based on 
#' proportion.
#' @param genes subset of the genes to use in PCA (`getTopHVGs`, etc)
#' @param gamma param to SCimplify
#' @param k.knn param to SCimplify
#' @param mode average or sum
#' @param ID SingleCellExperiment rowData column
#' @param verbose control logging
#' @importFrom SummarizedExperiment rowData `rowData<-`
#' @export
#' 
superCellMat <- function(sce, genes=NULL, prop=0.2, pca=TRUE, gamma=10, k.knn=5, rank=10,
    mode="average", ID="ID", verbose=TRUE) {
    if (!requireNamespace("SuperCell", quietly = TRUE)) {
        stop("SuperCell package needed to compute the metacell abundance. Please visit https://github.com/GfellerLab/SuperCell")
    } else {
        requireNamespace("SuperCell")
    }
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Matrix package needed to compute the metacell abundance. Please visit https://github.com/GfellerLab/SuperCell")
    } else {
        requireNamespace("Matrix")
        t <- Matrix::t
    }
    if (is.null(genes)) {
        genes <- getTopHVGs(sce, prop=prop)
    }

    GE <- sce@assays@data$logcounts
    row.names(GE) <- rowData(sce)[[ID]]

    if (pca) {
        pcs <- stats::prcomp(t(GE[genes, ]), rank. = rank)$x
        SC <- SuperCell::SCimplify_from_embedding(
          X = pcs, 
          k.knn = k.knn, 
          gamma = gamma
        )
    } else {
        SC <- SuperCell::SCimplify_from_embedding(
          X = GE, 
          k.knn = k.knn, 
          gamma = gamma
        )        
    }
    SC.GE <- SuperCell::supercell_GE(GE, SC$membership, mode=mode)
    cat("  ", dim(SC.GE),"\n")
    ret = SingleCellExperiment::SingleCellExperiment(assays=S4Vectors::SimpleList("logcounts"=SC.GE))
    rowData(ret) <- rowData(sce)
    return(ret)
}