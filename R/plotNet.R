#' plotNet
#' plot inferred network (with coefficients if available)
#' @param net network in bn class
#' @export
plotNet <- function(net, data=NULL, layout="kk", geom=geom_edge_link,
	largest_components=TRUE, sizeRange=c(1,5)) {
	  edgeArgs <- list()
	  edgeArgs[["color"]] <- "grey80"
	  edgeArgs[["arrow"]] <- arrow(length=unit(2,"mm"), type="closed")
	  edgeArgs[["start_cap"]] <- circle(5,"mm")
	  edgeArgs[["end_cap"]] <- circle(5,"mm")
  

	  geom <- geom_edge_link

	  if (is.null(data)) {
			g <- net %>% as.igraph() %>% as_tbl_graph()
	  } else {
	  	g <- bn.fit(net, data) %>% bn_fit_to_igraph() %>% as_tbl_graph()
	  	edgeArgs[["mapping"]] <- aes(color=coef)
	  	  		  edgeArgs[["color"]] <- NULL

	  }

if (largest_components) {
	comps <- to_components(g)
	g <- comps[[which.max(lapply(comps, function(x) {d <- data.frame(x) %>% dim(); d[1]}))]]
}

gra <- g %>%  mutate(degree=centrality_degree(mode="all")) %>%
ggraph(layout=layout) +
  do.call(geom, edgeArgs) +
  geom_node_point(aes(size=degree, color=degree)) +
  geom_node_text(aes(label=name, size=degree), repel=TRUE, bg.colour="white")+
  scale_size(range=sizeRange)+
  theme_graph()

if (is.null(data)) {
	gra
} else {
	gra + scale_edge_color_gradient(low="blue", high="red")
}

}