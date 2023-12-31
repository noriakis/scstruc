#' globalStrucValues
globalStrucValues <- function(spe, global_tbl_graph=NULL, labels, exclude_label=NA,
	summarize_func=mean, variation_func=sd, bn=NULL) {
    if (is.null(bn) & is.null(global_tbl_graph)) {stop("Please supply either of tbl_graph or bn")}
    if (!is.null(bn)) {return(.globalStrucValuesCoef(spe, bn, labels, exclude_label))}
    gedges <- global_tbl_graph |> activate(edges) |> data.frame()
    gnnames <- global_tbl_graph |> activate(nodes) |> pull(name)
    gedges$from <- gnnames[gedges$from]
    gedges$to <- gnnames[gedges$to]
    meta <- colData(spe) |> data.frame()    
    
    alllabels <- unique(colData(spe)[[label]])
    alllabels <- alllabels[!alllabels %in% exclude_label]

    logc <- spe@assays@data$logcounts
    if ("Symbol" %in% colnames(rowData(spe))) {
        row.names(logc) <- rowData(spe)$Symbol
    }
    ## If bn is not specified, the sum of abundances of two nodes connected by the edges are returned.
    ## Currently, multiple labels are not supported for this operation.
    if (length(labels)>1) {stop("Multiple labels are supported for coefficient mode")}
    appendix <- do.call(cbind, lapply(alllabels, function(x) {
        if ("barcode_id" %in% (meta |> colnames())) {
            inc_cells <- meta[meta[[label]] == x, ]$barcode_id
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

.globalStrucValuesCoef <- function(spe, bn, labels, exclude_label) {
    logc <- spe@assays@data$logcounts
    ## In case
    if ("Symbol" %in% colnames(rowData(spe))) {
        row.names(logc) <- rowData(spe)$Symbol
    }
    meta <- colData(spe) |> data.frame()
	if (!"barcode_id" %in% (meta |> colnames())) {
	    meta[["barcode_id"]] <- row.names(meta)
	}
	
	get_edges <- function(df) {## df means group.by metadata
	    group_name <- df[,labels] |> as.matrix()
	    group_name <- group_name[1,] |> as.character()
	    fit_df <- logc[names(bn$nodes),
	      which(colData(spe)$Barcode %in% df$Barcode)] |>
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