#' superCellMat
#' 
#' Return the metacell abundance matrix
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
superCellMat <- function(GE, pca=FALSE, genes=NULL, gamma=10, k.knn=5, mode="average") {
    if (pca) {
        GE <- stats::prcomp(Matrix::t(GE[genes, ]), rank. = 10)$x
        SC <- SCimplify_from_embedding(
          X = GE, 
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
    return(SC.GE)
}