#' perLabelStruc
#' 
#' Infer the BN per label (like using only the Epithelial cells or one sample)
#' The input is SingleCellExperiment or SpatialExperiment object and candidate genes to be included in the inference,
#' and label column in colData.
#' 
#' @param spe SingleCellExperiment or SpatialExperiment
#' @param all infer all the network of the corresponding label
#' @param label_name if all=FALSE, must be specified
#' @export
perLabelStruc <- function(spe, label, candidate_genes, algorithm="mmhc", reg=FALSE, penalty="glmnet",
    use_assay="logcounts", all=FALSE, label_name=NULL, verbose=FALSE, algorithm.args=list(), barcode_column="Barcode",
    change_symbol=TRUE, symbol_column="Symbol", cluster_label=NULL) {
    if (all) {
        alll <- unique(colData(spe)[[label]])
        nets <- lapply(alll, function(x) {
            cat(x, "\n")
            input <- .getInput(spe, candidate_genes, label, x, use_assay, barcode_column, cluster_label, verbose, change_symbol, symbol_column)
            net <- .getStruc(input, algorithm, reg, penalty, algorithm.args, verbose)
            return(net)
        })
        names(nets) <- alll
        return(nets)        
    } else {
        if (is.null(label_name)) {stop("Please specify label name")}
        input <- .getInput(spe, candidate_genes, label, label_name, use_assay, barcode_column, cluster_label, verbose,
            change_symbol, symbol_column)
        net <- .getStruc(input, algorithm, reg, penalty, algorithm.args, verbose)
        return(net)
    }
}

#' @noRd
.getInput <- function(sce, candidate_genes, label, x, use_assay, barcode_column="Barcode", cluster_label=NULL, verbose=FALSE,
    change_symbol=TRUE, symbol_column="Symbol") {
    logc <- sce@assays@data[[use_assay]]
    if (change_symbol) {
        row.names(logc) <- rowData(sce)[[symbol_column]]
    }
    meta <- colData(sce) |> data.frame()
    if (barcode_column=="row") {
        meta[["Barcode"]] <- row.names(meta)
        barcode_column <- "Barcode"
    }
    if (is.null(colnames(logc))) {
        colnames(logc) <- meta[[barcode_column]]
    }
    if (is.null(x)) {
        inc_cells <- meta[[barcode_column]]
    } else {
        inc_cells <- meta[meta[[label]] == x, ][[barcode_column]]
    }
    if (!is.null(cluster_label)) {
        ## First subset to inc_cell, and subsequently aggregate based on `cluster_label`.
        sce <- sce[, inc_cells %in% colData(sce)[[barcode_column]]]
        agg <- aggregateAcrossCells(sce, colData(sce)[[cluster_label]], use.assay.type=use_assay)
        input <- agg@assays@data[[use_assay]][intersect(row.names(agg), candidate_genes), ] |>
        as.matrix() |> t() |> 
        data.frame(check.names=FALSE)
    } else {
        input <- logc[intersect(row.names(logc), candidate_genes),
                  intersect(colnames(logc), inc_cells)] |>
        as.matrix() |> t() |>
        data.frame(check.names=FALSE)        
    }
    if (verbose) {
        cat("Dimension of the input for structure learning: n", dim(input)[1], "p", dim(input)[2], "\n")
    }
    if (dim(input)[2]==0) {stop("No genes")}
    if (dim(input)[1]==0) {stop("No samples")}
    return(input)
}


.getStruc <- function(input, algorithm, reg, penalty, algorithm.args, verbose) {
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
            algorithm.args[["data"]] <- input
            # algorithm.args[["restrict"]] <- "mmpc"
            # algorithm.args[["maximize"]] <- "hc"
            algorithm.args[["penalty"]] <- penalty
            net <- do.call(skeleton.reg, algorithm.args)
            return(net)
        }
    } else {
        algorithm.args[["x"]] <- input
        net <- do.call(algorithm, algorithm.args)
        return(net)
    }
}