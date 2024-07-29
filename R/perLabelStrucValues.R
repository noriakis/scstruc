# #' perLabelStrucValues
# #' 
# #' @param strucs named list of bn.strength
# #' 
# #' @export
# #' 
# perLabelStrucValues <- function(strucs, remove_zero=TRUE) {
#     booted_sub <- lapply(names(strucs), function(b) {
#         str <- strucs[[b]] %>% data.frame()
#         str[["coefficient"]] <- str[["strength"]]
#         str <- str[,c("from","to","coefficient")]
#         str[["label"]] <- b
#         str <- str[str[["coefficient"]]>0, ]
#         str
#     })
#     names(booted_sub) <- names(strucs)
#     merge_str <- do.call(rbind, booted_sub)
#     row.names(merge_str) <- 1:nrow(merge_str)
#     return(merge_str)
# }

#' @noRd
# DEPRECATED_perLabelStrucValues <- function(strucs) {
# 	logc <- spe@assays@data$logcounts
# 	alllabels <- unique(colData(spe)[[label]])
#     meta <- colData(spe) |> data.frame()

# 	converted <- lapply(seq_len(length(alllabels)), function(x) {

# 	    if (!is.null(per_struc[[x]])) {
# 	        graph <- per_struc[[x]] |> bnlearn::as.igraph() |> tidygraph::as_tbl_graph() |>
# 	            activate("edges") |> mutate(group=alllabels[x])      
# 	        g2 <- graph |> activate(edges) |> data.frame()
# 	        nnames <- graph |> activate(nodes) |> pull(name)
# 	        g2$from <- nnames[g2$from]
# 	        g2$to <- nnames[g2$to]
	        
# 	        inc_cells <- meta[meta[[label]] == alllabels[x], ]$barcode_id
# 	        frommat <- logc[g2$from, intersect(colnames(logc), inc_cells)]
# 	        tomat <- logc[g2$to, intersect(colnames(logc), inc_cells)]
# 	        avemat <- apply(frommat + tomat, 1, mean)
# 	        sdmat <- apply(frommat + tomat, 1, sd)
	        
# 	        g2$edgeval <- avemat |> as.numeric()
# 	        g2$edgesd <- sdmat |> as.numeric()
	        
# 	        tbl_graph(edges=g2)
# 	    }
#     })

#     comparegraph <- Reduce(graph_join,
#        converted[lapply(converted, function(x) !is.null(x)) |> unlist()])
#     return(comparegraph)
# }