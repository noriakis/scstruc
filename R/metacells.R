#' superCellMat
#' 
#' Return the metacell abundance matrix
#' Should provide SingleCellExperiment with logcounts data filled.
#' 
#' @import SuperCell
#' @param GE gene (row) to cell (column) matrix, like log-normalized
#' @param pca use PCA or not
#' @param genes subset of the genes to use in PCA (`getTopHVGs`, etc)
#' @param gamma param to SCimplify
#' @param k.knn param to SCimplify
#' @param mode average or sum
#' @export
#' 
superCellMat <- function(sce, genes=NULL, prop=0.2, pca=TRUE, gamma=10, k.knn=5, rank=10,
    mode="average", ID="ID", verbose=TRUE) {
    if (is.null(genes)) {
        genes <- getTopHVGs(sce, prop=prop)
    }

    GE <- sce@assays@data$logcounts
    row.names(GE) <- rowData(sce)[[ID]]

    if (pca) {
        pcs <- stats::prcomp(Matrix::t(GE[genes, ]), rank. = rank)$x
        SC <- SCimplify_from_embedding(
          X = pcs, 
          k.knn = k.knn, 
          gamma = gamma
        )
    } else {
        SC <- SCimplify_from_embedding(
          X = GE, 
          k.knn = k.knn, 
          gamma = gamma
        )        
    }
    SC.GE <- supercell_GE(GE, SC$membership, mode=mode)
    cat("  ", dim(SC.GE),"\n")
    ret = SingleCellExperiment(assays=SimpleList("logcounts"=SC.GE))
    rowData(ret) <- rowData(sce)
    return(ret)
}