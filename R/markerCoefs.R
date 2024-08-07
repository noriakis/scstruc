#' @title markerCoefs
#' @description identify marker regulatory relationships based on 
#' fitted coefficient matrix.
#' @details Users can choose XGBoost or Boruta algorithm.
#' @param coef_mat returned value of `strucValues`
#' @param classif_label which label to be used for classification
#' @param cell_label cell label for fitting
#' @param cell_column which column represents cell types
markerCoefs <- function(coef_mat, classif_label="group",
	cell_label=NULL, cell_column="label", sample_column="Sample",
    tentative_fix=TRUE, return_mat=FALSE, verbose=FALSE, returnChar=TRUE,
    xgboost=FALSE, xgboostArgs=list()) {
	if (is.null(cell_label)) {
		cell_label <- coef_mat[[cell_column]] |> unique()
	}
    mydata <- coef_mat |>
        dplyr::mutate(edge_name=paste0(from,"->",to)) |>
        dplyr::filter(.data[[cell_column]] %in% cell_label) |>
        tidyr::pivot_wider(id_cols=edge_name,
            values_from=coefficient,
            names_from=sample_column) |>
        data.frame()
    row.names(mydata) <- mydata$edge_name
    mydata$edge_name <- NULL
    mydata <- mydata |> t() |> data.frame(check.names=FALSE)

    group_dic <- coef_mat[,c(sample_column, classif_label)]
    group_dic <- group_dic[!duplicated(group_dic),]

    group_dic <- group_dic[[classif_label]] |> setNames(group_dic[[sample_column]])

    mydata[["classif_label"]] <- group_dic[row.names(mydata)] |> factor()

    if (return_mat) {
        return(mydata)
    }

    ## We do not use impXgboost in Boruta as not recommended.
    if (xgboost) {
        require(xgboost)
        trueLabel <- unique(as.character(mydata$classif_label))[2]
        cat("True label will be:", trueLabel, "\n")
        vec <- as.numeric(mydata$classif_label==trueLabel)
        mydata$classif_label <- NULL
        if (length(xgboostArgs)==0) {
            xgboostArgs[["nrounds"]] <- 100
        }
        xgboostArgs[["data"]] <- mydata %>% as.matrix()
        xgboostArgs[["label"]] <- vec
        res <- do.call(xgboost::xgboost, xgboostArgs)
        importance <- xgb.importance(feature_names = colnames(mydata), model = res)
        return(list("xgboost"=res, "importance"=importance, "data"=mydata))
    } else {
        mydata[is.na(mydata)] <- 0
        cat("Performing Boruta algorithm ...\n")
        brt <- Boruta::Boruta(classif_label ~ ., data=mydata)
        if (tentative_fix) {
            brt_fixed <- TentativeRoughFix(brt)
        } else {
            brt_fixed <- brt
        }
        important_edge <- brt_fixed$finalDecision[
            brt_fixed$finalDecision=="Confirmed"
        ] |>
            names()
        if (!returnChar) {
            important_edge <- do.call(rbind, important_edge %>% strsplit("->")) %>% data.frame() %>% `colnames<-`(c("from","to"))
        }
        return(list(brt_fixed, important_edge))    
    }
}
