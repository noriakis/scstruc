#' globalStrucValues
globalStrucValues <- function(spe, global_tbl_graph, label, exclude_label=NA,
	summarize_func=mean, variation_func=sd) {
    gedges <- global_tbl_graph |> activate(edges) |> data.frame()
    gnnames <- global_tbl_graph |> activate(nodes) |> pull(name)
    gedges$from <- gnnames[gedges$from]
    gedges$to <- gnnames[gedges$to]
    meta <- colData(spe) |> data.frame()    
    
    alllabels <- unique(colData(spe)[[label]])
    alllabels <- alllabels[!alllabels %in% exclude_label]
    
    logc <- spe@assays@data$logcounts
    
    appendix <- do.call(cbind, lapply(alllabels, function(x) {
        inc_cells <- meta[meta[[label]] == x, ]$barcode_id
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
