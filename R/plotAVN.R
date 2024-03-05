#' plotAVN
#' @param ew edge width, strength or direction
#' @export
plotAVN <- function(avn, str, edge=geom_edge_diagonal, ew=strength, lyt="fr",
    size_degree=TRUE) {
  
  ## Find undirected arcs
  avdf <- igraph::as_data_frame(bnlearn::as.igraph(avn))
  g <- igraph::graph_from_data_frame(merge(avdf, str, by=c("from","to")))

  # if (dim(ud)[1]!=0) {
  #   g <- g %>% as_tbl_graph() %E>% mutate(und=from %in% ud[,1] & to %in% ud[,2])
  # } else {
  #   g <- g %>% as_tbl_graph() %E>% mutate(und=FALSE)
  # }  
  V(g)$degree <- igraph::degree(g)
  
  gg <- ggraph(g, layout=lyt)
  # gg <- gg + do.call(edge, list("mapping"=aes(width=!!enquo(ew), color=!!enquo(ew), filter=und)))
  ## Directed arcs
  gg <- gg + do.call(edge, list("mapping"=aes(width=!!enquo(ew), color=!!enquo(ew)),
                                "arrow"=arrow(length=unit(2,"mm"), type="closed"),
                                "start_cap"=circle(5,"mm"),
                                "end_cap"=circle(5,"mm")))
  if (size_degree) {
      gg <- gg + geom_node_point(aes(size=degree, color=degree))
  } else {
      gg <- gg + geom_node_point()
  }
  if (size_degree) {
      gg <- gg + geom_node_text(aes(label=name, size=degree, color=degree), repel=TRUE, bg.colour="white")  
  } else {
      gg <- gg + geom_node_text(aes(label=name), repel=TRUE, bg.colour="white")
  }
  gg <- gg + theme_graph()
  gg <- gg + scale_edge_width(range=c(0.1, 2))
  gg <- gg + scale_edge_color_gradient(low="blue",high="red")
  gg <- gg + scale_color_gradient(low="blue",high="red")
  return(gg)
}
