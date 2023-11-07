markerCoefs <- function(coef_mat, spe, classif_label="group",
	cell_label=NULL, cell_column="label", sample_column="Sample.Name"
    tentative_fix=TRUE) {
	if (is.null(cell_label)) {
		cell_label <- coef_mat[[cell_column]] |> unique()
	}
    mydata <- coef_mat |>
        mutate(edge_name=paste0(from,"->",to)) |>
        filter(.data[[cell_column]] %in% cell_label) |>
        tidyr::pivot_wider(id_cols=edge_name,
            values_from=coefficient,
            names_from=sample_column) |>
        data.frame()
    row.names(mydata) <- mydata$edge_name
    mydata$edge_name <- NULL
    mydata <- mydata |> t() |> data.frame(check.names=FALSE)
    library(Boruta)
    brt <- Boruta(classif_label ~ ., data=mydata)
    if (tentative_fix) {
        brt_fixed <- TentativeRoughFix(brt)
    } else {
        brt_fixed <- brt
    }
    important_edge <- brt_fixed$finalDecision[
        brt_fixed$finalDecision=="Confirmed"
    ] |>
        names()
    return(list(brt_fixed, important_edge))
}