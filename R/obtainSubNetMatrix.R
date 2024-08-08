#' obtainSubNetMatrix
#' 
#' @param df the results of `strucValues` (specifying the bn object)
obtainSubNetMatrix <- function(df, spe, func="sum", use_assay="logcounts") {
	logc <- spe@assays@data[[use_assay]]

	all_nodes <- unique(c(df$from, df$to))
	mat <- do.call(rbind, lapply(all_nodes, function(x) {
	    subset_graph <- df |> filter(
	        from %in% x | to %in% x
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