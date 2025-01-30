#' @title plotNetPerCell
#' @description The function is for plotting the coefficient inside the graph per cell.
#' @param values the data.frame output of strucValues(), storing per-cell (per-cluster) parameters.
#' If needed, values should be subset beforehand for the interesting clusters.
#' @param net igraph or bn object
#' @param avn averaged network
#' @param str strength object
#' @param lyt layout
#' @param size_degree whether to size the nodes using the centrality_degree()
#' @export
plotNetPerCell <- function(
    values,
    net=NULL,
    avn=NULL,
    str=NULL,
    lyt="fr",
    size_degree=TRUE
    ) {
    if (!is.null(avn) & !is.null(str)) {
        ## Find undirected arcs
        avdf <- igraph::as_data_frame(bnlearn::as.igraph(avn))
        g <- igraph::graph_from_data_frame(merge(avdf, str, by=c("from","to")))
        igraph::V(g)$degree <- igraph::degree(g)
    } else {
        if (is.null(net)) {
            stop("Please provide net (igraph or bn object), or avn (averaged network of bnlearn) and str (strength from bnlearn)")
        }
        if ("igraph" %in% class(net)) {
            g <- net   
        } else {
            g <- bnlearn::as.igraph(net)
        }
        igraph::V(g)$degree <- igraph::degree(g)
    }
    g2 <- as_tbl_graph(g)
    g2 <- g2 %>% bind_edges(values)
    ## Merge the values per cell
    
  
    gg <- ggraph(g2, layout=lyt)
    ew <- "coefficient"
    # gg <- gg + do.call(edge, list("mapping"=aes(width=!!enquo(ew), color=!!enquo(ew), filter=und)))
    ## Directed arcs
    gg <- gg + geom_edge_parallel(mapping=aes(width=.data$coefficient,
    	color=.data$label, filter=!is.na(.data$coefficient)),
        arrow=arrow(length=unit(2,"mm"), type="closed"),
        start_cap=circle(5,"mm"),
        end_cap=circle(5,"mm"))
    if (size_degree) {
        gg <- gg + geom_node_point(aes(size=.data$degree, color=.data$degree))
    } else {
        gg <- gg + geom_node_point()
    }
    if (size_degree) {
        gg <- gg + geom_node_text(aes(label=.data$name,
        	size=.data$degree, color=.data$degree), repel=TRUE, bg.colour="white")  
    } else {
        gg <- gg + geom_node_text(aes(label=.data$name), repel=TRUE, bg.colour="white")
    }
    gg <- gg + theme_graph()
    gg <- gg + scale_edge_width(range=c(0.1, 5))
    # gg <- gg + scale_edge_color_gradient(low="blue",high="red")
    # gg <- gg + scale_color_gradient(low="blue",high="red")
    return(gg)
}
