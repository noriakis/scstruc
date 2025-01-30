#' @title obtainSubNetMatrix
#' @description obtain subnetwork abundance matrix
#' @param df the results of `strucValues` (specifying the bn object)
#' @param spe SingleCellExperiment object
#' @param func how to summarize subnetwork abundance
#' @param useAssay which assay to use
#' @export
obtainSubNetMatrix <- function(df, spe, func="sum", useAssay="logcounts") {
	logc <- spe@assays@data[[useAssay]]

	all_nodes <- unique(c(df$from, df$to))
	mat <- do.call(rbind, lapply(all_nodes, function(x) {
	    subset_graph <- df |> filter(
	        .data$from %in% x | .data$to %in% x
	    )
	    subset_nodes <- unique(c(subset_graph$from, subset_graph$to))
	    in_net <- logc[subset_nodes, ]
	    if (is.vector(in_net)) {
	    	return(in_net)
	    } else {
    	    return(apply(in_net, 2, func))
	    }
	}))
	row.names(mat) <- paste0(all_nodes, "_net")
	return(mat)
}