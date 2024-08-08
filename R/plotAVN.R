#' @title plotAVN
#' @description plot the averaged network based on bootstrapping
#' @param str boot.strength object
#' @param edge geom to be used in edge
#' @param ew edge width, strength or direction
#' @param layout layout in ggraph
#' @param sizeDegree size the node based on degree
#' @param sizeRange node size range
#' @param degreeMode degree mode in igraph::degree
#' @param largestComponents take largest components automatically
#' @param threshold threshold for averaged network
#' @export
plotAVN <- function(str, edge=geom_edge_link,
    ew=strength, layout="kk", degreeMode="all",
    sizeDegree=TRUE, sizeRange=c(3, 6),
    largestComponents=TRUE, threshold=NULL) {
    
    if (is.null(threshold)) {
        avn <- bnlearn::averaged.network(str)
    } else {
        avn <- bnlearn::averaged.network(str, threshold=threshold)
    }
    ## Find undirected arcs
    avdf <- igraph::as_data_frame(bnlearn::as.igraph(avn))
    g <- igraph::graph_from_data_frame(merge(avdf, str, by=c("from","to")))

    V(g)$degree <- igraph::degree(g, mode=degreeMode)

    g <- g %>% as_tbl_graph()

    if (largestComponents) {
        comps <- to_components(g)
        g <- comps[[which.max(lapply(comps, function(x) {d <- data.frame(x) %>% dim(); d[1]}))]]
    }

    gg <- ggraph(g, layout=layout)
    ## Directed arcs
    gg <- gg + do.call(edge, list("mapping"=aes(width=!!enquo(ew), color=!!enquo(ew)),
                                "arrow"=arrow(length=unit(2,"mm"), type="closed"),
                                "start_cap"=circle(5,"mm"),
                                "end_cap"=circle(5,"mm")))
    if (sizeDegree) {
        gg <- gg + geom_node_point(aes(size=degree, color=degree))
    } else {
        gg <- gg + geom_node_point()
    }
    if (sizeDegree) {
        gg <- gg + geom_node_text(aes(label=name, size=degree), repel=TRUE, bg.colour="white")  
    } else {
        gg <- gg + geom_node_text(aes(label=name), repel=TRUE, bg.colour="white")
    }
    gg <- gg + theme_graph()
    gg <- gg + scale_size(range=sizeRange)
    gg <- gg + scale_edge_width(range=c(0.1, 2))
    gg <- gg + scale_edge_color_gradient()
    gg <- gg + scale_color_gradient()
    return(gg)
}
