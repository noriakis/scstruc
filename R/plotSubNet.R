#' @title plotSubNet
#' @description plot subnetwork centered to specified gene
#' 
#' @param df the results of `strucValues`
#' @param candidate_node_id candidate node id to subset
#' @param cell_label if you need to plot per cell basis, specify cell label,
#' along with cell column
#' @param cell_column if you need to plot per cell basis, specify cell column,
#' along with cell label
#' @param coef_cutoff edge with value below this thrshold will not be drawn.
#' @param with_hist histology image plotting
#' @param spe SpatialExperiment class object
#' @param label name in spe
#' @param node_size node size
#' @param show which label to show
#' @param layout layout in ggraph
#' @param cell_label cell label name in cell_column
#' @param cell_column cell column name
#' @param coef_cutoff coefficient cutoff to be plotted
#' @param sort sort the label
#' @param highlight_edge edge highlighting
#' @return ggplot
#' @export
plotSubNet <- function(df, candidate_node_id, sort=FALSE,
    with_hist=FALSE, spe=NULL, label="label", node_size=3, show=NULL,
    layout="kk", cell_label=NULL, cell_column=NULL, highlight_edge=NULL,
    coef_cutoff=0.01) {
    # Probably better to split to sce and spatial experiment
    flag <- is.null(cell_label)+is.null(cell_column)
    if (flag==1) {stop("Please specify both of cell_label and cell_column")}

    if (with_hist) {
        if (is.null(spe)) {stop("Please provide SpatialExperiment object")}
        myst <- plotST(spe, label=label)
        cols <- ggplot_build(myst)$data[[2]] |> dplyr::select("fill", "group")
        cols <- cols[!duplicated(cols),]
        alllabels <- factor(myst$data[[4]]) |> levels()
        colvec <- cols$fill |> setNames(alllabels[cols$group])
    }
    if (is.null(show)) {show <- df[[label]] |> unique()}
    subset_graph <- df |> 
        filter(.data[[label]] %in% show)

    if (!is.null(cell_label)) {
        cat("Subsetting to ", cell_label, "\n")
        subset_graph <- subset_graph |> 
            filter(.data[[cell_column]] %in% cell_label)
    }

    subset_graph <- subset_graph |>  
        filter(
            .data$from %in% candidate_node_id | .data$to %in% candidate_node_id
        )
    if (dim(subset_graph)[1]==0) {stop("No edges present")}
    minc <- subset_graph$coefficient |> min(na.rm=TRUE)
    maxc <- subset_graph$coefficient |> max(na.rm=TRUE)
    subset_graph$edge_name <- paste0(subset_graph$from,"->",subset_graph$to)

    plot_subset_g <- tbl_graph(edges=subset_graph)

    labels <- plot_subset_g |>
        activate("edges") |>
        pull(!!enquo(label)) |>
        unique()
    if (sort) {labels <- labels |> sort()}

    plot_list <- lapply(labels,
        function (x) {
            if (with_hist) {
                curcol <- colvec[x]           
            } else {
                curcol <- "black"
            }
            returnplot <- plot_subset_g |>
            activate("edges") |>
            filter(.data[[label]] == !!enquo(x)) |>
            activate("nodes") |>
            mutate(group = x) |>
            ggraph(layout=layout) +
            geom_node_point(color=curcol, size=node_size) +
            geom_edge_parallel(
                aes(color=.data$coefficient, width=.data$coefficient, filter= abs(.data$coefficient) > coef_cutoff),
                end_cap=circle(2,"mm"),
                start_cap=circle(2,"mm"),
                arrow=arrow(length=unit(1.5,"mm"), type="closed"))+
            ggfx::with_outer_glow(
            geom_edge_parallel(
                aes(color=.data$coefficient, width=.data$coefficient,
                filter=.data$edge_name %in% highlight_edge & abs(.data$coefficient) > coef_cutoff),
                end_cap=circle(2,"mm"),
                start_cap=circle(2,"mm"),
                arrow=arrow(length=unit(1.5,"mm"), type="closed")),
                colour="gold", sigma=5
            )+
            scale_edge_color_gradient2(low="blue", mid="white", high="red", limits=c(minc, maxc))+
            scale_edge_width(range=c(0.5, 1.5), limits=c(minc, maxc))+
            geom_node_text(aes(label=.data$name), repel=TRUE, bg.colour="white")+
            ggtitle(paste0(x))+ theme_graph()
            if (!is.null(cell_label)) {
                returnplot <- returnplot + labs(subtitle=paste0("(", cell_label,")"))
            }
            returnplot
        })
    names(plot_list) <- labels
    wrapped <- patchwork::wrap_plots(plot_list, guides="collect")
    if (with_hist) {
        patchwork::wrap_plots(wrapped, myst+theme_graph())
    } else {
        wrapped
    }
}
