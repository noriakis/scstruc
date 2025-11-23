
#' When the excessive zero, some folds contain only the value zero
#' which throws an error in cv.ncvreg and cv.glmnet.
#' This function ignores the error and returns no candidate when the
#' function raises an error.
#' @noRd
.ncvregCV <- function(X, y, nFolds, algorithm, s=NULL) {
    if (!requireNamespace("ncvreg")) {
        stop("Needs ncvreg installation")
    }
    algorithm <- strsplit(algorithm, "_") %>% vapply("[", 1, FUN.VALUE="a")
    tryCatch({
        fit <- ncvreg::cv.ncvreg(X, y, alpha=1, penalty=algorithm, nfolds=nFolds)
        included_var_index <- as.numeric(which(fit$fit$beta[, fit$min][-1]!=0))
        included_var_index
        return(included_var_index)
        },
        error = function(e) {
            message(e)
            message("If constant error, return no candidate for this node\n")                
            return(NULL)
        },
        finally={}
    )

}

#' @noRd
#' @importFrom glmnet cv.glmnet
.glmnetCV <- function(X, y, nFolds=5, algorithm=NULL,
    s="lambda.min", standardize=TRUE) {

    tryCatch({
        fit <- cv.glmnet(X %>% as.matrix(), y,
            alpha=1, family="gaussian", nfolds=nFolds,
            standardize=standardize)
        numcoef <- coef(fit, s=s)[,1]
        numcoef <- numcoef[2:length(numcoef)]
        included_var_index <- which(numcoef!=0)
        return(included_var_index)
        },
        error = function(e) {
            message(e)
            message("If constant error, return no candidate for this node\n")                
            return(NULL)
        },
        finally={}
    )
}


#' @noRd
#' @importFrom glmnet glmnet
.glmnetBIC <- function(X, y, nFolds=NULL, algorithm=NULL,
    s=NULL, onlyBIC=FALSE) {
    ## This corresponds to L1MB in Schmidt et al. 2007. AAAI.
    fit <- glmnet(X, y, alpha=1, family="gaussian")
    tLL <- fit$nulldev - fit$nulldev * (1-fit$dev.ratio)
    k <- fit$df
    n <- nobs(fit)
    BIC <- log(n)*k - tLL
    if (onlyBIC) {
        bicdf <- data.frame(lambda=fit$lambda, BIC=BIC)
        bicplot <- ggplot(bicdf, aes(x=.data$lambda, y=.data$BIC)) + 
            geom_point() +  
            geom_hline(color="red", yintercept=BIC[which.min(BIC)])+
            cowplot::theme_cowplot()
        return(
            list("fit"=fit,
                "data"=bicdf,
                "plot"=bicplot)
        )
    }
    included_var_index <- which(as.numeric(fit$beta[,which.min(BIC)])!=0)
    included_var_index
}

#' glmnetBICpath
#' return glmnet BIC path for lambdas
#' @export
#' @param data data for fitting
#' @param nn node name for Y
#' @examples
#' library(ggplot2)
#' data(gaussian.test)
#' glmnetBICpath(gaussian.test, "A")$plot
glmnetBICpath <- function(data, nn) {
    nodes <- colnames(data)
    X <- data[, setdiff(nodes, nn)] %>% as.matrix()
    y <- data[, nn]
    .glmnetBIC(X, y, onlyBIC=TRUE)
}

# .glmmTMB <- function(nodes, nn, data, formula= ~ 1) {
#     pred <- setdiff(nodes, nn)
#     if (length(pred) == 0) {
#         model = paste(nn, "~ 1")
#     } else {
#         model = paste(nn, "~", paste(pred, collapse = "+"))
#     }
#     fit <- glmmTMB(as.formula(model), data=data,
#         ziformula= formula,
#         family=gaussian)
#     res <- summary(fit)
#     tbl <- data.frame(res$coefficients$cond)
#     tbl <- tbl[2:nrow(tbl),]
#     tbl <- tbl[!is.na(tbl[,4]),]
#     row.names(tbl[tbl[,4]<0.05,])
# }

#' @noRd
.plasso <- function(X, y, lambda=0.1, nFolds=NULL, algorithm=NULL, s=NULL) {
    fit <- plasso_fit(X %>% as.matrix(), y, lambda=lambda, maxIter=100, gamma=0.5, eps=1e-6)
    ## The last coefficient is intercept
    included_var_index <- which(fit[1:(length(fit)-1)]!=0)
    included_var_index
}


#' For reproducibility, CV folds are created by caret
#' (Not in use)
#' @noRd
# returnFolds <- function(y, nFolds) {
#     flds <- caret::createFolds(y,
#         k=nFolds, list=TRUE, returnTrain=FALSE)
#     foldids = rep(0,length(y))
#     for (i in seq_len(nFolds)) {
#         unval <- unique(y[flds[[paste0("Fold",i)]]])
#         if (length(unval) == 1) {message("A fold containing constants!")}
#         foldids[flds[[paste0("Fold",i)]]] = i
#     }
#     return(foldids)
# }


#' @noRd
.L0CV <- function(X, y, nFolds=5, algorithm=NULL, s=NULL) {
    if (!requireNamespace("L0Learn")) {
        stop("Needs installation of L0Learn")
    } else {
        requireNamespace("L0Learn")
    }
    tryCatch({
        fit=L0Learn::L0Learn.cvfit(x=X %>% as.matrix(), y=y, penalty="L0",
            nFolds=nFolds, algorithm="CD")
        lambdaIndex <- which.min(fit$cvMeans[[1]])
        cand_lambda <- fit$fit$lambda[[1]][lambdaIndex]
        L0coef <- coef(fit$fit, lambda=cand_lambda)
        vars <- as.numeric(L0coef)[2:length(L0coef)]
        included_var_index <- which(vars!=0)
        return(included_var_index)
    }, error=function(e) {message(e); return(NULL)})
}

#' @noRd
.L0LXCV <- function(X, y, nFolds=5, algorithm=NULL, s=NULL) {
    if (!requireNamespace("L0Learn")) {
        stop("Needs installation of L0Learn")
    } else {
        requireNamespace("L0Learn")
    }
    algorithm <- strsplit(algorithm, "_") %>% vapply("[", 1, FUN.VALUE="a")
    tryCatch({
        fit=L0Learn::L0Learn.cvfit(x=X, y=y, penalty=algorithm, nFolds=nFolds, algorithm="CD")
        ## Choose gamma which have the lowest PE
        min_cvmean <- which.min(lapply(fit$cvMeans, function(cvm) cvm[[which.min(cvm)]]) |> unlist())
        chosen_gamma <- fit$fit$gamma[[min_cvmean]]

        ## Lambda in the selected gamma
        lambdaList <- fit$fit$lambda[[min_cvmean]]

        ## the lambda with the lowest PE
        chosen_lambda <- lambdaList[which.min(fit$cvMeans[[min_cvmean]])]
        L0coef <- coef(fit$fit, gamma=chosen_gamma, lambda=chosen_lambda)[,1]
        vars <- as.numeric(L0coef)[2:length(L0coef)]
        included_var_index <- which(vars!=0)
        return(included_var_index)        
    }, error=function(e) {message(e); return(NULL)})

}

#' @noRd
.MASTHurdle <- function(nodes, nn, data, formula= ~ 1) {
    pred <- setdiff(nodes, nn)
    if (length(pred) == 0) {
        model = paste(nn, "~ 1")
    } else {
        model = paste(nn, "~", paste(pred, collapse = "+"))
    }
    fit <- MAST::zlm(as.formula(model), sca=data)
    if (is.null(fit$cont)) {return(NULL)}
    res <- summary.glm(fit$cont)
    tbl <- data.frame(coefficients(res))
    tbl <- tbl[2:nrow(tbl),]
    tbl <- tbl[!is.na(tbl[,4]),]
    return(row.names(tbl[tbl[,4]<0.05,]))
}

#' @title skeleton.reg.boot
#' @description bootstrapping the two-stage approaches
#' @param data input data
#' @param algorithm two-stage algorithm
#' @param R replicate number
#' @param m sampling number
#' @param algorithm.args passed to algorithm function
#' @export
skeleton.reg.boot <- function(data, algorithm="glmnet_CV",
    R=100,  m=nrow(data), algorithm.args=list()) {
    nodes = names(data)
    perRun <- list()
    for (r in seq_len(R)) {
        resampling = sample(nrow(data), m, replace = TRUE)
        # generate the r-th bootstrap sample.
        replicate = data[resampling, , drop = FALSE]
        algorithm.args[["data"]] <- replicate
        algorithm.args[["algorithm"]] <- algorithm
        run <- do.call(skeleton.reg, algorithm.args)
        perRun[[r]] <- run
    }    
    st <- custom.strength(perRun, nodes)
    return(st)
}

#' Custom score function for bnlearn, using hurdle model in glmmTMB.
#' BIC is scaled by -2.
# glmmtmb.bic <- function(node, parents, data, args) {
#     ## Using glmmTMB
#     if (length(parents) == 0) {
#         model = paste(node, "~ 1")
#     } else {
#         model = paste(node, "~", paste(parents, collapse = "+"))
#     }
#     # cat("Fitting:", model, "\n")

#     fit <- glmmTMB(as.formula(model), data=data,
#         ziformula= ~ 1, family=gaussian, control = glmmTMBControl(parallel = 1))

#    - BIC(fit) / 2

# }#MY.BIC

#' Custom score function for bnlearn, using hurdle model.
#' @importFrom MAST zlm
#' @param node candidate node
#' @param parents parents of the nodes
#' @param data data used for scoring
#' @param args other arguments
#' @export
hurdle.bic <- function(node, parents, data, args) {
    if (isTRUE(args$cdrAdjustment)) {
        cdr <- rowSums(data>0)
        data[["cdr"]] <- scale(cdr)
    }
    if (length(parents)==0) {
        if (isTRUE(args$cdrAdjustment)) {
            mod <- paste0(node, " ~ cdr")
        } else {
            mod <- paste0(node, " ~ 1")
        }
    } else {
        if (isTRUE(args$cdrAdjustment)) {
            mod <- paste0(node, "~", paste0(parents, collapse=" + "), " + cdr")
        } else {
            mod <- paste0(node, "~", paste0(parents, collapse=" + "))
        }
    }
    # cat("Fitting ZLM:", model, "\n")
    fit <- MAST::zlm(as.formula(mod), sca=data)
    if (is.null(fit$cont)) return(-Inf)
    bic.sum <- -1 * (BIC.zlm.bayesglm(llk.zlm.bayesglm(fit$cont)) + 
        BIC.zlm.bayesglm(llk.zlm.bayesglm(fit$disc)))
    return(bic.sum)
}



#' @noRd
llk.zlm.bayesglm <- function(obj) {
  fam <- obj$family$family
  p <- obj$rank
  if(fam %in% c("gaussian", "Gamma", "inverse.gaussian")) p <- p + 1
  val <- p - obj$aic / 2
  attr(val, "nobs") <- sum(!is.na(obj$residuals))
  attr(val, "df") <- p
  class(val) <- "logLik"
  val
}


#' @noRd
BIC.zlm.bayesglm <- function(obj) {
  -2 * as.numeric(obj) + attr(obj, "df") * log(nobs(obj))
}


#' Custom score function for bnlearn, using hurdle model.
#' @noRd
#' @importFrom MAST zlm
hurdle.aic <- function(node, parents, data, args) {

    if (args$cdrAdjustment) {
        cdr <- rowSums(data>0)
        data[["cdr"]] <- scale(cdr)
    }
    if (length(parents)==0) {
        if (args$cdrAdjustment) {
            mod <- paste0(node, "~ cdr")
        } else {
            mod <- paste0(node, " ~ 1")
        }
    } else {
        if (args$cdrAdjustment) {
            mod <- paste0(node, "~", paste0(parents, collapse=" + "), " + cdr")
        } else {
            mod <- paste0(node, "~", paste0(parents, collapse=" + "))
        }
    }

    fit <- MAST::zlm(as.formula(mod), sca=data)
    if (is.null(fit$cont)) return(-Inf)    
    aic.sum <- -1 * (fit$cont$aic + fit$disc$aic)
    return(aic.sum)

}