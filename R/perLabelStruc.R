perLabelStruc <- function(spe, label, candidate_genes, algorithm="mmhc", reg=FALSE) {
	logc <- spe@assays@data$logcounts
	meta <- colData(spe) |> data.frame()
    nets <- lapply(unique(colData(spe)[[label]]), function(x) {
	    inc_cells <- meta[meta[[label]] == x, ]$barcode_id
	    input <- logc[candidate_genes,
	                  intersect(colnames(logc), inc_cells)] |>
	        as.matrix() |> t() |>
	        data.frame(check.names=FALSE)
	    
	    if (dim(input)[2]==0) {return(NULL)}
	    if (dim(input)[1]==0) {return(NULL)}
	    if (reg) {
		    pergraph <- bnlearnReg::rsmax2(input, restrict="mmpc", maximize="hc",
		              penalty="glmnet", nFolds=5, debug=F)
	    } else {
	        pergraph <- do.call(algorithm, list(x=input))	
	    }
	    return(pergraph)
	})
    return(nets)
}