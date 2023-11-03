#' perLabelStrucValues
perLabelStrucValues <- function(spe, per_struc, label) {
	logc <- spe@assays@data$logcounts
	alllabels <- unique(colData(spe)[[label]])
    meta <- colData(spe) |> data.frame()

	converted <- lapply(seq_len(length(alllabels)), function(x) {

	    if (!is.null(per_struc[[x]])) {
	        graph <- per_struc[[x]] |> bnlearn::as.igraph() |> tidygraph::as_tbl_graph() |>
	            activate("edges") |> mutate(group=alllabels[x])      
	        g2 <- graph |> activate(edges) |> data.frame()
	        nnames <- graph |> activate(nodes) |> pull(name)
	        g2$from <- nnames[g2$from]
	        g2$to <- nnames[g2$to]
	        
	        inc_cells <- meta[meta[[label]] == alllabels[x], ]$barcode_id
	        frommat <- logc[g2$from, intersect(colnames(logc), inc_cells)]
	        tomat <- logc[g2$to, intersect(colnames(logc), inc_cells)]
	        avemat <- apply(frommat + tomat, 1, mean)
	        sdmat <- apply(frommat + tomat, 1, sd)
	        
	        g2$edgeval <- avemat |> as.numeric()
	        g2$edgesd <- sdmat |> as.numeric()
	        
	        tbl_graph(edges=g2)
	    }
    })

    comparegraph <- Reduce(graph_join,
       converted[lapply(converted, function(x) !is.null(x)) |> unlist()])
    return(comparegraph)
}