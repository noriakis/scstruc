#' @noRd
.getInput <- function(sce,
    candidate_genes, label, x,
    use_assay,
    cluster_label=NULL, verbose=FALSE,
    change_symbol=TRUE, symbol_column="Symbol",
    zeroFilt=1, rmNeg=FALSE) {
    logc <- sce@assays@data[[use_assay]]
    if (change_symbol) {
        row.names(logc) <- rowData(sce)[[symbol_column]]
    }
    # if (is.null(colData(sce)[[barcode_column]])) {
    #     colData(sce)[[barcode_column]] <- as.character(1:ncol(sce))
    # }

    meta <- colData(sce) |> data.frame()

    # if (barcode_column=="row") {
    #     if (is.null(row.names(colData(sce)))) {
    #          row.names(meta) <- as.character(1:ncol(sce))
    #     }
    # }

    # if (barcode_column=="row") {
    #     meta[["Barcode"]] <- row.names(meta)
    #     barcode_column <- "Barcode"
    # }
    # if (is.null(colnames(logc))) {
    # colnames(logc) <- meta[[barcode_column]]
    # }
    if (is.null(x)) {
        inc_cells <- seq_len(ncol(sce))
    } else {
        inc_cells <- which(meta[[label]] %in% x)
    }
    if (!is.null(cluster_label)) {
        ## First subset to inc_cell, and subsequently aggregate based on `cluster_label`.
        sce <- sce[, inc_cells]# %in% colData(sce)[[barcode_column]]]
        agg <- aggregateAcrossCells(sce, colData(sce)[[cluster_label]], use.assay.type=use_assay)
        input <- agg@assays@data[[use_assay]][intersect(row.names(agg), candidate_genes), ] |>
        as.matrix() |> t() |> 
        data.frame(check.names=FALSE)
    } else {
        input <- logc[intersect(row.names(logc), candidate_genes),
                  inc_cells] |>
        as.matrix() |> t() |>
        data.frame(check.names=FALSE)        
    }
    if (dim(input)[1]==1) {stop("Only one gene is available")}
    if (verbose) {
        cat_subtle("Dimension of the input: n ", dim(input)[1], ", p ", dim(input)[2], "\n")
    }
    if (dim(input)[2]==0) {stop("No genes remained")}
    if (dim(input)[1]==0) {stop("No samples remained")}
    input <- input[, apply(input==0, 2, function(x) sum(x) <= nrow(input) * zeroFilt)]
    if (rmNeg) {
        input <- removeAllNegative(input)
    }
    # input <- input[, apply(input, 2, function(x) unique(x)!=1)]
    if (verbose) {
        cat_subtle("Dimension of the input for structure learning (filtered by the number of zero): n ",
            dim(input)[1], ", p ", dim(input)[2], "\n")
    }
    return(input)
}

