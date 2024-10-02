#' sid.bn
#' @description given the bn object, calculate SID between them using the original implementation.
#' @param true.bn true bn object
#' @param est.bn estimated bn object
#' @param twoway return the sid(true.bn, est.bn) + sid(est.bn, true.bn) divided by 2
#' @export
sid.bn <- function(true.bn, est.bn, twoway=FALSE) {
  if (!requireNamespace("SID")) {
    stop("Needs installation of SID")
  } else {
    requireNamespace("SID")
  }
  if (twoway) {
      res1 <- SID::structIntervDist(as_adjacency_matrix(as.igraph(true.bn)),
                               as_adjacency_matrix(as.igraph(est.bn)))
      res2 <- SID::structIntervDist(as_adjacency_matrix(as.igraph(est.bn)),
                               as_adjacency_matrix(as.igraph(true.bn)))
      return((res1$sid+res2$sid)/2)

  } else {
      res <- SID::structIntervDist(as_adjacency_matrix(as.igraph(true.bn)),
                                   as_adjacency_matrix(as.igraph(est.bn)))
      return(res$sid)

  }
}