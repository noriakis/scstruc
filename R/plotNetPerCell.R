#' plotNetPerCell
#' @param values the data.frame output of globalStrucValues(), storing per-cell (per-cluster) parameters.
#' If needed, values should be subset beforehand for the interesting clusters.
#' @param net igraph or bn object
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
        V(g)$degree <- igraph::degree(g)
    } else {
        if (is.null(net)) {
            stop("Please provide net (igraph or bn object), or avn (averaged network of bnlearn) and str (strength from bnlearn)")
        }
        if ("igraph" %in% class(net)) {
            g <- net        
        } else {
            g <- bnlearn::as.igraph(net)
        }
        V(g)$degree <- igraph::degree(g)
    }
    g2 <- as_tbl_graph(g)
    g2 <- g2 %>% bind_edges(values)
    ## Merge the values per cell
    
  
    gg <- ggraph(g2, layout=lyt)
    ew <- "coefficient"
    # gg <- gg + do.call(edge, list("mapping"=aes(width=!!enquo(ew), color=!!enquo(ew), filter=und)))
    ## Directed arcs
    gg <- gg + geom_edge_parallel(mapping=aes(width=coefficient, color=label, filter=!is.na(coefficient)),
        arrow=arrow(length=unit(2,"mm"), type="closed"),
        start_cap=circle(5,"mm"),
        end_cap=circle(5,"mm"))
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
    gg <- gg + scale_edge_width(range=c(0.1, 5))
    # gg <- gg + scale_edge_color_gradient(low="blue",high="red")
    # gg <- gg + scale_color_gradient(low="blue",high="red")
    return(gg)
}
