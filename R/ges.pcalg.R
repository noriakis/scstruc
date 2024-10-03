#' ges.pcalg
#' Perform Greedy Equivalent Search implemented in pcalg
#' @param data data
#' @param args arguments to pcalg::ges
#' @export
#' @importFrom pcalg ges
ges.pcalg <- function(data, args=NULL) {
	if (is.null(args)) {
		args <- list()
	}
	score <- new(pcalg::.__C__GaussL0penObsScore, data=data)
	args[["score"]] <- score
	args[["iterate"]] <- TRUE
	ges.fit <- do.call(pcalg::ges, args)

	g <- ges.fit$repr$weight.mat()
	ig <- igraph::graph_from_adjacency_matrix(g, mode="directed",weighted = TRUE, diag=TRUE)
	bn <- bnlearn::as.bn(ig)
	return(bn)	
}