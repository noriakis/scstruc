#' @title plotNet
#' @description plot inferred network (with coefficients if available)
#' The nodes will be sized and colored based on degree calculated by
#' tidygraph::centrality_degree.
#' @param net network in bn class
#' @param data data to fit the inferred BN
#' @param layout layout of the graph
#' @param geom geom for edge in ggraph
#' @param largeset_components take largest components automatically
#' @param sizeRange specifying size of the nodes
#' @param edgeArgs edge arguments used in `geom` argument
#' @param degreeMode degree method
#' @param highlightEdges if the two column matrix (corresponding to from and to)
#' is provided, highlight the corresponding edges in the plot.
#' @param highlightEdgeWidth highlight edge width
#' @export
plotNet <- function(net, data=NULL, layout="kk", geom=geom_edge_link,
    largest_components=TRUE, sizeRange=c(3,6), edgeArgs=list(),
    degreeMode="all", highlightEdges=NULL, highlightEdgeWidth=1) {
      
    if (is.null(edgeArgs[["color"]])) {
          edgeArgs[["color"]] <- "grey80"
    }
    if (is.null(edgeArgs[["arrow"]])) {
          edgeArgs[["arrow"]] <- arrow(length=unit(2,"mm"), type="closed")
    }
    if (is.null(edgeArgs[["start_cap"]])) {
        edgeArgs[["start_cap"]] <- circle(5,"mm")
    }
    if (is.null(edgeArgs[["end_cap"]])) {
          edgeArgs[["end_cap"]] <- circle(5,"mm")
    }
  
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

    if (is.null(highlightEdges)) {
        gra <- g %>%  mutate(degree=centrality_degree(mode=degreeMode)) %>%
        ggraph(layout=layout) +
            do.call(geom, edgeArgs) +
            geom_node_point(aes(size=degree, color=degree)) +
            geom_node_text(aes(label=name, size=degree), repel=TRUE, bg.colour="white")+
            scale_size(range=sizeRange)+
            theme_graph()      
    } else {
        if (!is.matrix(highlightEdges)) {
            stop("Please provide two column matrix")
        }
        edgeArgs2 <- edgeArgs
        edgeArgs[["mapping"]] <- c(edgeArgs[["mapping"]], aes(filter=!highlight))
        edgeArgs2[["mapping"]] <- c(aes(filter=highlight))
        edgeArgs2[["color"]] <- "tomato"
        edgeArgs2[["width"]] <- highlightEdgeWidth
        nl <- g %N>% pull(name)
        g <- g %E>% mutate(fromn = nl[from], ton = nl[to]) %>%
          mutate(highlight=(fromn %in% highlightEdges[,1]) & (ton %in% highlightEdges[,2]))
        gra <- g %N>%  mutate(degree=centrality_degree(mode=degreeMode)) %>%
        ggraph(layout=layout) +
          do.call(geom, edgeArgs) +
          do.call(geom, edgeArgs2) +
          geom_node_point(aes(size=degree, color=degree)) +
          geom_node_text(aes(label=name, size=degree), repel=TRUE, bg.colour="white")+
          scale_size(range=sizeRange)+
          theme_graph()
    }


    if (is.null(data)) {
          return(gra)
    } else {
          return(gra + scale_edge_color_gradient())
    }

}