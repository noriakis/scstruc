markerNetwork <- function(calcmat, spe, label="label", positive_label="1",
    tentative_fix=TRUE) {
    if (!requireNamespace("Boruta")) {
        stop("Needs Boruta installation")
    }

    X <- calcmat |> t() |> data.frame(check.names = FALSE)
    if ("barcode_id" %in% colnames(colData(spe))) {
        Y <- colData(spe)[[label]] |>
            setNames(colData(spe)$barcode_id)        
    } else {
        Y <- colData(spe)[[label]] |>
            setNames(colData(spe) |> row.names())
    }
    Y <- Y[row.names(X)]
    X[["resp"]] <- ifelse(Y==positive_label, 1, 0)
    X <- X[!is.na(X[["resp"]]),]

    brt <- Boruta::Boruta(resp ~ ., data=X)
    if (tentative_fix) {
        brt_fixed <- Boruta::TentativeRoughFix(brt)
    } else {
        brt_fixed <- brt
    }
    important_net <- brt_fixed$finalDecision[
        brt_fixed$finalDecision=="Confirmed"
    ] |>
        names()
    return(list(brt_fixed, important_net))
}