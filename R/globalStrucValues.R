#' @title strucValues
#' @description Fit the network parameters using SingleCellExperiment object and bn object.
#' 
#' @export
strucValues <- function(spe, global_tbl_graph=NULL, labels, exclude_label=NA,
	summarize_func=mean, variation_func=sd, bn=NULL, barcode="barcode_id", assay="logcounts",
    verbose=FALSE) {
    if (is.null(bn) & is.null(global_tbl_graph)) {stop("Please supply either of tbl_graph or bn")}
    if (!is.null(bn)) {return(.strucValuesCoef(spe, bn, labels, exclude_label, barcode, assay, verbose))}
    gedges <- global_tbl_graph |> activate(edges) |> data.frame()
    gnnames <- global_tbl_graph |> activate(nodes) |> pull(name)
    gedges$from <- gnnames[gedges$from]
    gedges$to <- gnnames[gedges$to]
    meta <- colData(spe) |> data.frame()    
    
    alllabels <- unique(colData(spe)[[label]])
    alllabels <- alllabels[!alllabels %in% exclude_label]

    logc <- spe@assays@data[[assay]]
    if ("Symbol" %in% colnames(rowData(spe))) {
        row.names(logc) <- rowData(spe)$Symbol
    }
    ## If bn is not specified, the sum of abundances of two nodes connected by the edges are returned.
    ## Currently, multiple labels are not supported for this operation.
    if (length(labels)>1) {stop("Multiple labels are supported for coefficient mode")}
    appendix <- do.call(cbind, lapply(alllabels, function(x) {
        if (barcode %in% (meta |> colnames())) {
            inc_cells <- meta[meta[[label]] == x, ][[barcode]]
        } else {
            inc_cells <- meta[meta[[label]] == x, ] |> row.names()
        }
        frommat <- logc[gedges$from, intersect(colnames(logc), inc_cells)]
        tomat <- logc[gedges$to, intersect(colnames(logc), inc_cells)]
        avemat <- apply(frommat + tomat, 1, summarize_func)
        sdmat <- apply(frommat + tomat, 1, variation_func)
        cbind(avemat, sdmat) |>
            `colnames<-`(c(paste0(x,"_summarize"), paste0(x,"_var")))
    }))
    appendix <- cbind(gedges, appendix)
    row.names(appendix) <- seq_len(nrow(appendix))
    return(appendix)
}

#' @noRd
.strucValuesCoef <- function(spe, bn, labels, exclude_label, barcode, assay, verbose) {
    cat("Coefficient calculation per specified group:", paste0(labels, collapse=", "), "\n")
    logc <- spe@assays@data[[assay]]
    ## In case
    if ("Symbol" %in% colnames(rowData(spe))) {
        row.names(logc) <- rowData(spe)$Symbol
    }
    meta <- colData(spe) |> data.frame()
	if (!barcode %in% (meta |> colnames())) {
	    meta[[barcode]] <- row.names(meta)
        colData(spe)[[barcode]] <- row.names(meta)
	}

	get_edges <- function(df) {## df means group.by metadata

	    group_name <- df[,labels] |> as.matrix()
	    group_name <- group_name[1,] |> as.character()

	    fit_df <- logc[names(bn$nodes),
	      which(colData(spe)[[barcode]] %in% df[[barcode]])] |>
	      as.matrix() |> t() |>
	      data.frame(check.names=FALSE)

	    fitted <- bnlearn::bn.fit(bn, fit_df)
	    all_genes <- names(fitted)
	    all_genes_edges <- do.call(rbind, lapply(all_genes, function(gene) {
	        tmp_coef <- fitted[[gene]]$coefficients
	        coefn <- names(tmp_coef)
	        if (length(coefn)>1) {
	            coefn <- coefn[2:length(coefn)]
	            do.call(rbind, lapply(coefn, function(from) {
	                c(from, gene, tmp_coef[from], group_name)
	            }))
	        }
	    })) |> data.frame() |> `colnames<-`(c("from","to","coefficient", labels))
	    all_genes_edges$coefficient <- all_genes_edges$coefficient |> as.numeric()
	    return(all_genes_edges)
	}
	coefs <- meta |> group_by(across(all_of(labels))) |>
	    group_map(~get_edges(.x), .keep=TRUE)
    coefs <- do.call(rbind, coefs)
    return(coefs)
}