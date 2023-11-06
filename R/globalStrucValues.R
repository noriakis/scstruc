#' globalStrucValues
globalStrucValues <- function(spe, global_tbl_graph=NULL, label, exclude_label=NA,
    label2=NULL, exclude_label2=NA,
	summarize_func=mean, variation_func=sd, bn=NULL) {
    if (is.null(bn) & is.null(global_tbl_graph)) {stop("Please supply either of tbl_graph or bn")}
    if (!is.null(bn)) {return(.globalStrucValuesCoef(spe, bn, label, exclude_label, label2, exclude_label2))}
    gedges <- global_tbl_graph |> activate(edges) |> data.frame()
    gnnames <- global_tbl_graph |> activate(nodes) |> pull(name)
    gedges$from <- gnnames[gedges$from]
    gedges$to <- gnnames[gedges$to]
    meta <- colData(spe) |> data.frame()    
    
    alllabels <- unique(colData(spe)[[label]])
    alllabels <- alllabels[!alllabels %in% exclude_label]
    
    logc <- spe@assays@data$logcounts
    ## If bn is not specified, the sum of abundances of two nodes connected by the edges are returned.
    appendix <- do.call(cbind, lapply(alllabels, function(x) {
        if ("barcode_id" %in% (meta |> colnames())) {
            inc_cells <- meta[!meta[[label]] == x, ]$barcode_id
        } else {
            inc_cells <- meta[!meta[[label]] == x, ] |> row.names()
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

.globalStrucValuesCoef <- function(spe, bn, label, exclude_label, label2, exclude_label2) {
    
    alllabels <- unique(colData(spe)[[label]])
    alllabels <- alllabels[!alllabels %in% exclude_label]

    if (!is.null(label2)) {
        alllabels2 <- unique(colData(spe)[[label2]])
        alllabels2 <- alllabels2[!alllabels2 %in% exclude_label2]        
    }

    logc <- spe@assays@data$logcounts
    meta <- colData(spe) |> data.frame()
    if (length(intersect(names(bn$nodes), row.names(logc)))!=
        length(names(bn$nodes))) {
        stop("Gene name mismatch between bn and single-cell experiment")
    }

    coefs <- lapply(alllabels, function (x) {
        if ("barcode_id" %in% (meta |> colnames())) {
            inc_cells <- meta[!meta[[label]] == x, ]$barcode_id
        } else {
            inc_cells <- meta[!meta[[label]] == x, ] |> row.names()
        }
        if (is.null(label2)) {
            fit_df <- logc[names(bn$nodes),
              intersect(colnames(logc), inc_cells)] |>
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
                        c(from, gene, tmp_coef[from])
                    }))
                }
            })) |> data.frame() |> `colnames<-`(c("from","to","coefficient"))
            all_genes_edges$coefficient <- all_genes_edges$coefficient |> as.numeric()
            all_genes_edges$group <- x
            return(all_genes_edges)
        } else {
            ## If label2 is specified
            l2_edgelist <- do.call(rbind, lapply(alllabels2, function(l2) {
                if ("barcode_id" %in% (meta |> colnames())) {
                    inc_cells2 <- meta[!meta[[label2]] == l2, ]$barcode_id
                } else {
                    inc_cells2 <- meta[!meta[[label2]] == l2, ] |> row.names()
                }
                inc_cells_combined <- intersect(inc_cells2, inc_cells)
                fit_df <- logc[names(bn$nodes),
                    intersect(colnames(logc), inc_cells_combined)] |>
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
                            c(from, gene, tmp_coef[from])
                        }))
                    }
                })) |> data.frame() |> `colnames<-`(c("from","to","coefficient"))
                all_genes_edges$coefficient <- all_genes_edges$coefficient |> as.numeric()
                all_genes_edges$group2 <- l2
                return(all_genes_edges)
            }))
            l2_edgelist$group <- x
            return(l2_edgelist)
        }
        
    })
    names(coefs) <- alllabels
    coefs <- do.call(rbind, coefs)
    return(coefs)
}