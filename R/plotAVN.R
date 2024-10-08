#' @title plotAVN
#' @description plot the averaged network based on bootstrapping
#' @param str boot.strength object
#' @param edge geom to be used in edge
#' @param ew edge width, strength, direction, or coef
#' if coef, data to be fitted must be provided
#' @param ec edge color, strength, direction, or coef
#' if coef, data to be fitted must be provided
#' @param data data to fit the paramter on averaged network
#' @param layout layout in ggraph
#' @param sizeDegree size the node based on degree
#' @param sizeRange node size range
#' @param degreeMode degree mode in igraph::degree
#' @param largestComponents take largest components automatically
#' @param threshold threshold for averaged network
#' @param highlightEdges highlight edges, two-column matrix is expected
#' If highlight is specified, ggfx::with_outer_glow 
#' will be used to highlight the edge
#' @param highlightColor highlight color for ggfx
#' @param highlightExpand expand parameter for ggfx
#' @param highlightSigma sigma paramter for ggfx
#' @importFrom ggplot2 arrow
#' @import ggraph
#' @export
plotAVN <- function(str, edge=geom_edge_link,
    ew="strength", ec="direction", layout="kk", degreeMode="all",
    sizeDegree=TRUE, sizeRange=c(3, 6), data=NULL,
    highlightColor="tomato", highlightExpand=5, highlightSigma=2,
    highlightEdges=NULL, largestComponents=TRUE, threshold=NULL) {
    if (!is.null(highlightEdges)) {
        if (!requireNamespace("ggfx")) {
            stop("Please install ggfx")
        }
    }
    if (as.character(ew)=="coef" & is.null(data)) {
        stop("if specified coef, please provide data")
    }
    if (as.character(ec)=="coef" & is.null(data)) {
        stop("if specified coef, please provide data")
    }
    if (is.null(threshold)) {
        avn <- bnlearn::averaged.network(str)
    } else {
        avn <- bnlearn::averaged.network(str, threshold=threshold)
    }
    if (!is.null(data)) {
        fitted <- bn.fit(avn, data)
        avdf <- igraph::as_data_frame(bn_fit_to_igraph(fitted))
    } else {
        ## Find undirected arcs
        avdf <- igraph::as_data_frame(bnlearn::as.igraph(avn))        
    }
    g <- igraph::graph_from_data_frame(merge(avdf, str, by=c("from","to")))
    V(g)$degree <- igraph::degree(g, mode=degreeMode)

    g <- g %>% as_tbl_graph()

    if (largestComponents) {
        comps <- to_components(g)
        g <- comps[[which.max(lapply(comps, function(x) {d <- data.frame(x) %>% dim(); d[1]}))]]
    }
    if (ew=="direction") {
        g <- g %E>% mutate(direction1=direction)
        ew <- "direction1"
    }
    if (ec=="direction") {
        g <- g %E>% mutate(direction1=direction)
        ec <- "direction1"
    }

    if (!is.null(highlightEdges)) {
        highchar <- paste0(highlightEdges[,1], "->", highlightEdges[,2])
        nl <- g %N>% pull(name)
        g <- g %E>% mutate(fromn = nl[from], ton = nl[to]) %>%
          mutate(highlight=paste0(fromn, "->", ton) %in% highchar)
    }

    gg <- ggraph(g, layout=layout)
    ## Directed arcs
    if (!is.null(highlightEdges)) {
        gg <- gg + ggfx::with_outer_glow(do.call(edge, list("mapping"=aes("width"=.data[[ew]], "color"=.data[[ec]], "filter"=.data[["highlight"]]),
                                    "arrow"=arrow(length=unit(2,"mm"), type="closed"),
                                    "start_cap"=circle(5,"mm"),
                                    "end_cap"=circle(5,"mm"))), colour=highlightColor, expand=highlightExpand, sigma=highlightSigma)
        gg <- gg + do.call(edge, list("mapping"=aes("width"=.data[[ew]], "color"=.data[[ec]], "filter"=!.data[["highlight"]]),
                                    "arrow"=arrow(length=unit(2,"mm"), type="closed"),
                                    "start_cap"=circle(5,"mm"),
                                    "end_cap"=circle(5,"mm")))

    } else {
        gg <- gg + do.call(edge, list("mapping"=aes("width"=.data[[ew]], "color"=.data[[ec]]),
                                    "arrow"=arrow(length=unit(2,"mm"), type="closed"),
                                    "start_cap"=circle(5,"mm"),
                                    "end_cap"=circle(5,"mm")))
    }

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
    if (ec=="direction1") {ec <- "direction"}
    if (ew=="direction1") {ew <- "direction"}
    gg <- gg + scale_edge_width(range=c(0.1, 2), name=ew)
    gg <- gg + scale_edge_color_gradient(name=ec)
    gg <- gg + scale_color_gradient()
    return(gg)
}
