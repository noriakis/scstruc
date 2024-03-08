#' @param listOfNets list of bn.strength
#' @export
coreEdges <- function(listOfNets, label_name) {
    tblg <- tbl_graph(edges=do.call(rbind, lapply(seq_along(listOfNets), function(x) {
      avn <- data.frame(igraph::as_edgelist(bnlearn::as.igraph(averaged.network(listOfNets[[x]][[label_name]]))))
      colnames(avn) <- c("from","to")
      df <- merge(avn, x[[label_name]], by=c("from","to"))
      if (!is.null(names(listOfNets)) {
        df[["net"]] <- names(listOfNets)[x]
      }
      return(df)
    })))
    return(tblg)
}
