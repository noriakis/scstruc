
.ncvregCV <- function(X, y, nFolds, algorithm, s=NULL) {
    algorithm <- strsplit(algorithm, "_") %>% vapply("[", 1, FUN.VALUE="a")
    fit <- cv.ncvreg(X, y, alpha=1, penalty=algorithm, nfolds=nFolds)
    included_var_index <- as.numeric(which(fit$fit$beta[, fit$min][-1]!=0))
    included_var_index
}

.glmnetBIC <- function(X, y, nFolds=NULL, algorithm=NULL, s=NULL) {
    ## This corresponds to L1MB
    fit <- glmnet(X, y, alpha=1, family="gaussian")
    tLL <- fit$nulldev - fit$nulldev * (1-fit$dev.ratio)
    k <- fit$df
    n <- nobs(fit)
    BIC <- log(n)*k - tLL
    included_var_index <- which(as.numeric(fit$beta[,which.min(BIC)])!=0)
    included_var_index
}


.glmmTMB <- function(nodes, nn, data, formula= ~ 1) {
    pred <- setdiff(nodes, nn)
    if (length(pred) == 0) {
        model = paste(nn, "~ 1")
    } else {
        model = paste(nn, "~", paste(pred, collapse = "+"))
    }
    fit <- glmmTMB(as.formula(model), data=data,
        ziformula= formula,
        family=gaussian)
    res <- summary(fit)
    tbl <- data.frame(res$coefficients$cond)
    tbl <- tbl[2:nrow(tbl),]
    tbl <- tbl[!is.na(tbl[,4]),]
    row.names(tbl[tbl[,4]<0.05,])
}

.plasso <- function(X, y, lambda=0.1, nFolds=NULL, algorithm=NULL, s=NULL) {
    cat_subtle("Precision lasso lambda: ", lambda, "\n")
    fit <- plasso_fit(X %>% as.matrix(), y, lambda=lambda, maxIter=100, gamma=0.5, eps=1e-6)
    ## The last coefficient is intercept
    included_var_index <- which(fit[1:(length(fit)-1)]!=0)
    included_var_index
}

.glmnetCV <- function(X, y, nFolds=5, algorithm=NULL, s="lambda.min") {
    fit <- cv.glmnet(X %>% as.matrix(), y,
        alpha=1, family="gaussian", nfolds=nFolds)
    numcoef <- coef(fit, s=s)[,1]
    numcoef <- numcoef[2:length(numcoef)]
    included_var_index <- which(numcoef!=0)
    included_var_index
}

.L0CV <- function(X, y, nFolds=5, algorithm=NULL, s=NULL) {
    fit=L0Learn::L0Learn.cvfit(x=X %>% as.matrix(), y=y, penalty="L0",
        nFolds=nFolds, algorithm="CD")
    lambdaIndex <- which.min(fit$cvMeans[[1]])
    cand_lambda <- fit$fit$lambda[[1]][lambdaIndex]
    L0coef <- coef(fit$fit, lambda=cand_lambda)
    vars <- as.numeric(L0coef)[2:length(L0coef)]
    included_var_index <- which(vars!=0)
    return(included_var_index)
}


.L0LXCV <- function(X, y, nFolds=5, algorithm=NULL, s=NULL) {
    algorithm <- strsplit(algorithm, "_") %>% vapply("[", 1, FUN.VALUE="a")
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
    included_var_index
}

.MASTHurdle <- function(nodes, nn, data, formula= ~ 1) {

    pred <- setdiff(nodes, nn)
    if (length(pred) == 0) {
        model = paste(nn, "~ 1")
    } else {
        model = paste(nn, "~", paste(pred, collapse = "+"))
    }
    fit <- MAST::zlm(as.formula(model), data=data)
    res <- summary.glm(fit$cont)
    tbl <- data.frame(coefficients(res))
    tbl <- tbl[2:nrow(tbl),]
    tbl <- tbl[!is.na(tbl[,4]),]
    row.names(tbl[tbl[,4]<0.05,])
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
#' @noRd
glmmtmb.bic <- function(node, parents, data, args) {
    ## Using glmmTMB
    if (length(parents) == 0) {
        model = paste(node, "~ 1")
    } else {
        model = paste(node, "~", paste(parents, collapse = "+"))
    }
    # cat("Fitting:", model, "\n")

    fit <- glmmTMB(as.formula(model), data=data,
        ziformula= ~ 1, family=gaussian, control = glmmTMBControl(parallel = 1))

   - BIC(fit) / 2

}#MY.BIC

#' Custom score function for bnlearn, using hurdle model.
#' @noRd
#' @importFrom MAST zlm
hurdle.bic <- function(node, parents, data, args) {

    if (length(parents) == 0) {
        model = paste(node, "~ 1")
    } else {
        model = paste(node, "~", paste(parents, collapse = "+"))
    }
    # cat("Fitting ZLM:", model, "\n")
    fit <- MAST::zlm(as.formula(model), sca=data)
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

    if (length(parents) == 0) {
        model = paste(node, "~ 1")
    } else {
        model = paste(node, "~", paste(parents, collapse = "+"))
    }
    # cat("Fitting ZLM:", model, "\n")

    fit <- MAST::zlm(as.formula(model), sca=data)
    if (is.null(fit$cont)) return(-Inf)    
    aic.sum <- -1 * (fit$cont$aic + fit$disc$aic)
    return(aic.sum)    

}

#' skeleton.reg
#' 
#' Two-stage approach for regularization-based structure learning.
#' Currently, this function is experimental as 
#' it uses non-exported functions in `bnlearn`.
#' 
#' @param s effective only in glmnet, lambda.min or lambda.1se
#' @export
skeleton.reg <- function(data, algorithm="glmnet_CV", whitelist=NULL, blacklist=NULL,
    nFolds=5, verbose=FALSE, s="lambda.min", maximize="hc", maximize.args=list()) {
    if (verbose) {
        cat_subtle("Algorithm: ", algorithm, "\n")
        cat_subtle("Input for structure learning: n ", dim(data)[1], "p ", dim(data)[2], "\n")
    }
    nodes <- colnames(data)
    penFun <- dplyr::case_when(
        algorithm=="glmnet_BIC" ~ ".glmnetBIC",
        algorithm=="plasso" ~ ".plasso",
        algorithm=="glmnet_CV" ~ ".glmnetCV",
        algorithm=="MCP_CV" ~ ".ncvregCV",
        algorithm=="SCAD_CV" ~ ".ncvregCV",
        algorithm=="L0L2_CV" ~ ".L0LXCV",
        algorithm=="L0L1_CV" ~ ".L0LXCV",
        algorithm=="L0_CV" ~ ".L0CV"
    )

    mb <- sapply(nodes, function(nn) {
        if (verbose) {cat_subtle(" ", nn)}
        X <- data[, setdiff(nodes, nn)] %>% as.matrix()
        y <- data[, nn]
        if (length(unique(y))==1) {return(NULL)}
        if (algorithm=="glmmTMB") {
            sign <- .glmmTMB(nodes, nn, data)
            if (length(sign)==0) {return(NULL)}
            if (verbose) {cat_subtle(" candidate: ", length(sign), "\n")}
            return(sign)
        } else if (algorithm=="MAST") {
            sign <- .MASTHurdle(nodes, nn, data)
            if (length(sign)==0) {return(NULL)}
            if (verbose) {cat_subtle(" candidate: ", length(sign), "\n")}
            return(sign)
        } else {
            included_var_index <- do.call(penFun, list("X"=X, "y"=y,
                "nFolds"=nFolds, "algorithm"=algorithm, "s"=s))            
            if (verbose) {cat_subtle(" candidate: ", length(included_var_index), "\n")}
            if (length(included_var_index)==0) {return(NULL)}
            return(colnames(X)[included_var_index])
        }
    })

    names(mb) <- nodes

    ## Neighbors
    mb <- lapply(nodes, function(nn) {
        candidate.neighbours <- mb[[nn]]
        if (length(candidate.neighbours) == 0) {
            return(list(mb = character(0), nbr = character(0)))
        }
        return(list(mb = mb[[nn]], nbr = candidate.neighbours))
    })
    names(mb) <- nodes
    for (node in nodes)
      mb[[node]]$mb <- bnlearn:::fake.markov.blanket(mb, node)
    mb <- bnlearn:::bn.recovery(mb, nodes = nodes)

    arcs <- bnlearn:::nbr2arcs(mb)
    res <- list(learning = list("method"="scstruc",
        "blacklist"=list(blacklist)),
        nodes = bnlearn:::cache.structure(names(mb), arcs = arcs),
        arcs = arcs, blacklist=blacklist)
    constraints <- bnlearn:::arcs.to.be.added(res$arcs,
        nodes, whitelist = res$learning$blacklist)
    start <- empty.graph(nodes = nodes)
    maximize.args[["x"]] <- data
    maximize.args[["start"]] <- start
    maximize.args[["blacklist"]] <- constraints
    struc_res <- do.call(maximize, maximize.args)
    # hc_res <- hc(data, start=start, blacklist=constraints)
    return(struc_res)
}