#' plotNet
#' plot inferred network (with coefficients if available)
#' @param net network in bn class
#' @export
plotNet <- function(net, data=NULL, layout="kk", geom=geom_edge_link,
	largest_components=TRUE) {
	  edgeArgs <- list()
	  edgeArgs[["color"]] <- "grey80"
	  edgeArgs[["arrow"]] <- arrow(length=unit(2,"mm"), type="closed")
	  edgeArgs[["start_cap"]] <- circle(5,"mm")
	  edgeArgs[["end_cap"]] <- circle(5,"mm")
  

	  geom <- geom_edge_link
      g <- net %>% as.igraph() %>% as_tbl_graph()

if (largest_components) {
	comps <- to_components(g)
	g <- comps[[which.max(lapply(comps, function(x) {d <- data.frame(x) %>% dim(); d[1]}))]]
}

g %>%  mutate(degree=centrality_degree(mode="all")) %>%
ggraph(layout=layout) +
  do.call(geom, edgeArgs) +
  geom_node_point(aes(size=degree, color=degree)) +
  geom_node_text(aes(label=name, size=degree), repel=TRUE, bg.colour="white")+
  theme_graph()

}