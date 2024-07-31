## Merged to scstruc.R

# #' perLabelStruc
# #' 
# #' Infer the BN per label (like using only the Epithelial cells or one sample)
# #' The input is SingleCellExperiment or SpatialExperiment object and candidate genes to be included in the inference,
# #' and label column in colData.
# #' 
# #' @param spe SingleCellExperiment or SpatialExperiment
# #' @param all infer all the network of the corresponding label
# #' @param label_name if all=FALSE, must be specified
# #' @param return_data return the data used to infer the network
# #' @export
# perLabelStruc <- function(spe, candidate_genes, label="label", algorithm="mmhc",
#     use_assay="logcounts", all=FALSE, label_name=NULL, verbose=FALSE, algorithm.args=list(), barcode_column="Barcode",
#     change_symbol=TRUE, symbol_column="Symbol", cluster_label=NULL, nonzero=1, return_data=FALSE, rmNeg=FALSE) {
#     if (all) {
#         alll <- unique(colData(spe)[[label]])
#         nets <- lapply(alll, function(x) {
#             cat(x, "\n")
#             input <- .getInput(spe, candidate_genes, label, x, use_assay, barcode_column,
#                 cluster_label, verbose, change_symbol, symbol_column, nonzero, rmNeg   )
#             net <- .getStruc(input, algorithm, algorithm.args, verbose)
#             if (return_data) {
#                 return(list("net"=net, "data"=input))
#             } else {
#                 return(net)
#             }
#         })
#         names(nets) <- alll
#         return(nets)        
#     } else {
#         if (is.null(label_name)) {stop("Please specify label name")}
#         nets <- lapply(label_name, function(x) {
#             cat(x, "\n")
#             input <- .getInput(spe, candidate_genes, label, x, use_assay, barcode_column, cluster_label, verbose,
#                 change_symbol, symbol_column, nonzero, rmNeg)
#             net <- .getStruc(input, algorithm, algorithm.args, verbose)
#             if (return_data) {
#                 return(list("net"=net, "data"=input))
#             } else {
#                 return(net)
#             }
#         })
#         names(nets) <- label_name
#         return(nets)
#     }
# }

