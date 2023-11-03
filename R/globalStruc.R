#' globalStruc
globalStruc <- function(spe, candidate_genes, label, exclude_label=NA, algorithm="mmhc",
	reg=FALSE, algorithm.args=list()) {
	logc <- spe@assays@data$logcounts
	meta <- colData(spe) |> data.frame()
    inc_cells <- meta[!meta[[label]] %in% exclude_label, ]$barcode_id
    input <- logc[candidate_genes,
              intersect(colnames(logc), inc_cells)] |>
    as.matrix() |> t() |>
    data.frame(check.names=FALSE)
    if (reg) {
	    global_graph <- bnlearnReg::rsmax2(input, restrict="mmpc", maximize="hc",
	              penalty="glmnet", nFolds=5, debug=F)
    } else {
    	algorithm.args[["x"]] <- input
        global_graph <- do.call(algorithm, algorithm.args)	
    }
    global_tbl_graph <- global_graph |> 
        bnlearn::as.igraph() |>
        tidygraph::as_tbl_graph()
    return(global_tbl_graph)
}