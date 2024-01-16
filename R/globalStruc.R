#' globalStruc
#' @param algorithm inference algorithm of BN, ignored if reg=TRUE
#' @param cluster_labels (experimental) aggregate to speed up computation
#' @param penalty if {reg}=TRUE, the penalty used can be specified, default to {glmnet} (lasso).
#' If you want to use CCDr algorithm, specify `ccdr` or `ccdr.boot`.
#' If you use CCDr algorithm, the {bn} object is returned.
globalStruc <- function(spe, candidate_genes, label, exclude_label=NA, algorithm="mmhc",
	reg=FALSE, algorithm.args=list(), return_bn=FALSE, return_data=FALSE, debug=FALSE,
    cluster_label=NULL, penalty="glmnet",verbose=FALSE) {
	logc <- spe@assays@data$logcounts
	meta <- colData(spe) |> data.frame()
    if ("barcode_id" %in% (meta |> colnames())) {
        inc_cells <- meta[!meta[[label]] %in% exclude_label, ]$barcode_id
    } else {
        inc_cells <- meta[!meta[[label]] %in% exclude_label, ] |> row.names()
    }
    
    if (!is.null(cluster_label)) {
        agg <- aggregateAcrossCells(spe, cluster_label, use.assay.type="logcounts")
        input <- agg@assays@data$logcounts[candidate_genes, ] |>
        as.matrix() |> t() |> 
        data.frame(check.names=FALSE)
    } else {
        input <- logc[candidate_genes,
                  intersect(colnames(logc), inc_cells)] |>
        as.matrix() |> t() |>
        data.frame(check.names=FALSE)        
    }

    if (verbose) {
    	cat("Dimension of the input for structure learning: n", dim(input)[1], "p", dim(input)[2], "\n")
    }

    if (reg) {
    	if (verbose) {
	    	cat("Penalty: ", penalty, " selected\n")		
    	}
    	if (penalty=="ccdr") {
    		cat("Returning the result of ccdr.run\n")
	        dat <- sparsebnUtils::sparsebnData(input %>% as.matrix(), type = "continuous")
    		algorithm.args[["data"]] <- dat
    		return(do.call(ccdrAlgorithm::ccdr.run, algorithm.args))
    	} else if (penalty=="ccdr.boot") {
    		cat("Returning the list of bn.strength\n")
    		algorithm.args[["data"]] <- input
    		return(do.call(scstruc::ccdr.boot, algorithm.args))    		
    	} else {
		    global_graph <- bnlearnReg::rsmax2(input, restrict="mmpc", maximize="hc",
		              penalty=penalty, nFolds=5, debug=debug)    		
    	}
    } else {
    	algorithm.args[["x"]] <- input
    	algorithm.args[["debug"]] <- debug
        global_graph <- do.call(algorithm, algorithm.args)	
    }
    if (return_bn) {
        if (return_data) {
            return(list(global_graph, input))
        } else {
            return(global_graph)
        }
    }
    global_tbl_graph <- global_graph |> 
        bnlearn::as.igraph() |>
        tidygraph::as_tbl_graph()
    return(global_tbl_graph)
}