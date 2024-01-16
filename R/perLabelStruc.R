
#' @param all infer all the network of the corresponding label
#' @param label_name if all=FALSE, must be specified
#' @export
perLabelStruc <- function(spe, label, candidate_genes, algorithm="mmhc", reg=FALSE, penalty="glmnet",
    use_assay="logcounts", all=FALSE, label_name=NULL, verbose=FALSE, algorithm.args=list()) {
    cluster_label <- NULL
    if (all) {
        nets <- lapply(unique(colData(spe)[[label]]), function(x) {
            input <- .getInput(spe, candidate_genes, label, x, use_assay, cluster_label, verbose)
            net <- .getStruc(input, algorithm, reg, penalty, algorithm.args, verbose)
            return(net)
        })
        return(nets)        
    } else {
        if (is.null(label_name)) {stop("Please specify label name")}
        input <- .getInput(spe, candidate_genes, label, label_name, use_assay, cluster_label, verbose)
        net <- .getStruc(input, algorithm, reg, penalty, algorithm.args, verbose)
        return(net)
    }
}

#' @noRd
.getInput <- function(sce, candidate_genes, label, x, use_assay, cluster_label=NULL, verbose=FALSE) {
    logc <- sce@assays@data[[use_assay]]
    meta <- colData(sce) |> data.frame()
    if (is.null(x)) {
        inc_cells <- colnames(sce)
    } else {
        if ("barcode_id" %in% (meta |> colnames())) {
            inc_cells <- meta[meta[[label]] == x, ]$barcode_id
        } else {
            inc_cells <- meta[meta[[label]] == x, ] |> row.names()
        }
    }
    if (!is.null(cluster_label)) {
        ## Ignoring `x`
        agg <- aggregateAcrossCells(sce, colData(sce)[[cluster_label]], use.assay.type=use_assay)
        input <- agg@assays@data[[use_assay]][candidate_genes, ] |>
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
            algorithm.args[["x"]] <- input
            algorithm.args[["restrict"]] <- "mmpc"
            algorithm.args[["maximize"]] <- "hc"
            algorithm.args[["penalty"]] <- penalty
            net <- do.call(bnlearnReg::rsmax2, algorithm.args)
            return(net)
        }
    } else {
        algorithm.args[["x"]] <- input
        net <- do.call(algorithm, algorithm.args)
        return(net)
    }
}