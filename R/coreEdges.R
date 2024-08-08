#' @title coreBootEdges
#' @description concatenate bootstrapped Bayesian network
#' @param listOfNets list of bn.strength
#' @export
coreBootEdges <- function(listOfNets) {
    tblg <- tbl_graph(edges=do.call(rbind, lapply(seq_along(listOfNets), function(x) {
      avn <- data.frame(igraph::as_edgelist(bnlearn::as.igraph(averaged.network(listOfNets[[x]]$net))))
      colnames(avn) <- c("from","to")
      df <- merge(avn, listOfNets[[x]]$net, by=c("from","to"))
      if (!is.null(names(listOfNets))) {
        df[["net"]] <- names(listOfNets)[x]
      } else {
        df[["net"]] <- as.character(seq_along(listOfNets))
      }
      return(df)
    })))
    return(tblg)
}

#' @title coreEdges
#' @description concatenate Bayesian network
#' @param listOfNets list of bn
#' @export
coreEdges <- function(listOfNets) {
    Reduce(igraph::intersection, lapply(listOfNets, function(x) bnlearn::as.igraph(x))) %>%
    as_tbl_graph()
}
