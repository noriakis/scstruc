#' plotSubNet
#' @param df the results of `globalStrucValues`
#' @param candidate_node_id candidate node id to subset
#' @export
plotSubNet <- function(df, candidate_node_id, sort=FALSE,
    with_hist=FALSE, spe=NULL, label="label", node_size=3) {
    if (with_hist) {
        if (is.null(spe)) {stop("Please provide SpatialExperiment object")}
        myst <- plotST(spe, label=label)
        cols <- ggplot_build(myst)$data[[2]] |> dplyr::select(fill, group)
        cols <- cols[!duplicated(cols),]
        alllabels <- factor(myst$data[[4]]) |> levels()
        colvec <- cols$fill |> setNames(alllabels[cols$group])
    }

    subset_graph <- df |> filter(
        from %in% candidate_node_id | to %in% candidate_node_id
    )
    minc <- subset_graph$coefficient |> min()
    maxc <- subset_graph$coefficient |> max()

    plot_subset_g <- tbl_graph(edges=subset_graph)

    labels <- plot_subset_g |>
        activate(edges) |>
        pull(group) |>
        unique()
    if (sort) {labels <- labels |> sort()}
    plot_list <- lapply(labels,
        function (x) {
            if (with_hist) {
                curcol <- colvec[x]           
            } else {
                curcol <- "black"
            }
            plot_subset_g |>
            activate(edges) |>
            filter(.data$group == x) |>
            activate(nodes) |>
            mutate(group = x) |>
            ggraph(layout="kk") +
            geom_node_point(color=curcol, size=node_size) +
            geom_edge_link(
                aes(color=coefficient, width=coefficient),
                end_cap=circle(2,"mm"),
                start_cap=circle(2,"mm"),
                arrow=arrow(length=unit(1.5,"mm"), type="closed"))+
            scale_edge_color_gradient2(low="blue", high="red", limits=c(minc, maxc))+
            scale_edge_width(range=c(0.5, 1.5), limits=c(minc, maxc))+
            geom_node_text(aes(label=name), repel=TRUE, bg.colour="white")+
            ggtitle(x)+ theme_graph()
        })
    wrapped <- patchwork::wrap_plots(plot_list)
    if (with_hist) {
        patchwork::wrap_plots(wrapped, myst+theme_graph())
    } else {
        wrapped
    }
}
