#' @noRd
.getInput <- function(sce, candidate_genes, label, x, use_assay, barcode_column="Barcode", cluster_label=NULL, verbose=FALSE,
    change_symbol=TRUE, symbol_column="Symbol", nonzero=1) {
    logc <- sce@assays@data[[use_assay]]
    if (change_symbol) {
        row.names(logc) <- rowData(sce)[[symbol_column]]
    }
    if (is.null(colData(sce)[[barcode_column]])) {
        colData(sce)[[barcode_column]] <- as.character(1:ncol(sce))
    }
    meta <- colData(sce) |> data.frame()
    if (barcode_column=="row") {
        meta[["Barcode"]] <- row.names(meta)
        barcode_column <- "Barcode"
    }
    # if (is.null(colnames(logc))) {
    colnames(logc) <- meta[[barcode_column]]
    # }
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
        cat("Dimension of the input: n", dim(input)[1], "p", dim(input)[2], "\n")
    }
    if (dim(input)[2]==0) {stop("No genes")}
    if (dim(input)[1]==0) {stop("No samples")}
    input <- input[, apply(input==0, 2, function(x) sum(x) < nrow(input) * nonzero)]
    ## Remove constants
    # input <- input[, apply(input, 2, function(x) unique(x)!=1)]
    if (verbose) {
        cat("Dimension of the input for structure learning (filtered by the number of zero): n", dim(input)[1], "p", dim(input)[2], "\n")
    }
    return(input)
}


#' @noRd
.getStruc <- function(input, algorithm, algorithm.args, verbose) {
    if (verbose) {
        cat("Algorithm:", algorithm, "selected\n")        
    }
    if (algorithm=="ccdr") {
        cat("Returning the result of ccdr.run\n")
        dat <- sparsebnUtils::sparsebnData(input %>% as.matrix(), type = "continuous")
        algorithm.args[["data"]] <- dat
        return(do.call(ccdrAlgorithm::ccdr.run, algorithm.args))
    } else if (algorithm=="ccdr.boot") {
        cat("Returning the list of bn.strength\n")
        algorithm.args[["data"]] <- input
        return(do.call(scstruc::ccdr.boot, algorithm.args))         
    } else if (algorithm == "Hurdle") {
        algorithm.args[["data"]] <- input
        return(do.call(.Hurdle, algorithm.args))
    } else if (algorithm %in% c("mmhc","rsmax2","h2pc","hc","tabu")) {
        cat("Using default bnlearn algorithm\n")
        algorithm.args[["x"]] <- input
        net <- do.call(algorithm, algorithm.args)
        return(net)
    } else {
        if (grepl("boot", algorithm)) {
            cat("Bootstrapping specified\n")
            algorithm <- strsplit(algorithm, "\\.") %>% sapply("[", 1)
            nodes <- names(input)
            perRun <- list()
            m <- nrow(input) ## m determined
            if (is.null(algorithm.args[["R"]])) {
                cat("R not specified in bootstrapping, default to 100\n")
                algorithm.args[["R"]] <- 100
            }
            R <- algorithm.args[["R"]]
            algorithm.args[["R"]] <- NULL
            for (r in seq_len(R)) {
                if (verbose) {cat(r)}
                resampling = sample(nrow(input), m, replace = TRUE)
                replicate = input[resampling, , drop = FALSE]
                algorithm.args[["data"]] <- replicate
                algorithm.args[["algorithm"]] <- algorithm
                repnet <- do.call(skeleton.reg, algorithm.args)
                perRun[[r]] <- repnet
            }
            if (verbose) {cat("\n")}
            net <- custom.strength(perRun, nodes)                
        } else {
            algorithm.args[["verbose"]] <- verbose
            algorithm.args[["data"]] <- input
            # algorithm.args[["restrict"]] <- "mmpc"
            # algorithm.args[["maximize"]] <- "hc"
            algorithm.args[["algorithm"]] <- algorithm
            net <- do.call(skeleton.reg, algorithm.args)                
        }
        return(net)
    }
}