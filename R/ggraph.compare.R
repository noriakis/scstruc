#' ggraph.compare
#' compare two bn object based on parallel edges
#' @param nets (named) list of bn object
#' @param lyt layout in ggraph
#' @return ggplot2 object
#' @export
#' @examples
#' data(gaussian.test)
#' test <- head(gaussian.test, 10)
#' nets <- list("hc"=hc(test), "mmhc"=mmhc(test))
#' # ggraph.compare(nets)
ggraph.compare <- function(nets, lyt="kk") {
    if (is.null(names(nets))) {
        gn <- paste0("graph", seq_len(length(nets)))
        names(nets) <- gn
    } else {
        gn <- names(nets)
    }
    nets <- lapply(gn, function(x) {as.igraph(nets[[x]]) %>%
            as_tbl_graph() %E>% mutate(group=x)})
    concat <- Reduce(graph_join, nets)
    ggraph(concat, layout=lyt) +
        geom_edge_parallel(aes(color=.data$group),
                           arrow=arrow(type="closed", length = unit(1,"mm")),
                           start_cap=circle(2,"mm"), end_cap=circle(2,"mm"))+
        geom_node_point() +
        geom_node_text(aes(label=.data$name), repel=TRUE, bg.colour="white")+
        theme_graph()
}
