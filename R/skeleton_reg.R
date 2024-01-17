
#' @export
skeleton.reg.boot <- function(data, penalty="glmnet_CV", R=100,  m=nrow(data), algorithm.args=list()) {
    nodes = names(data)
    perRun <- list()
    for (r in seq_len(R)) {
        resampling = sample(nrow(data), m, replace = TRUE)
        # generate the r-th bootstrap sample.
        replicate = data[resampling, , drop = FALSE]
        algorithm.args[["data"]] <- replicate
        algorithm.args[["penalty"]] <- penalty
        run <- do.call(skeleton.reg, algorithm.args)
        perRun[[r]] <- run
    }    
    st <- custom.strength(perRun, nodes)
    return(st)
}




#' skeleton.reg
#' 
#' Two-stage approach for regularization-based structure learning.
#' Currently, this function is experimental as 
#' it uses non-exported functions in `bnlearn`.
#' 
#' @export
skeleton.reg <- function(data, penalty="glmnet_CV", whitelist=NULL, blacklist=NULL,
    nFolds=5) {
    nodes <- colnames(data)
    penFun <- dplyr::case_when(
        penalty=="glmnet_BIC" ~ ".glmnetBIC",
        penalty=="plasso_BIC" ~ ".plassoBIC",
        penalty=="glmnet_CV" ~ ".glmnetCV",
        penalty=="MCP_CV" ~ ".ncvregCV",
        penalty=="SCAD_CV" ~ ".ncvregCV",
        penalty=="L0L2_CV" ~ ".L0LXCV",
        penalty=="L0L1_CV" ~ ".L0LXCV",
        penalty=="L0_CV" ~ ".L0CV"
    )
    mb <- sapply(nodes, function(nn) {
        X <- data[, setdiff(nodes, nn)] %>% as.matrix()
        y <- data[, nn]
        included_var_index <- do.call(penFun, list("X"=X, "y"=y,
            "nFolds"=nFolds, "pen"=penalty))
        if (length(included_var_index)==0) {return(NULL)}
        colnames(X)[included_var_index]
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
    res <- list(learning = list("method"="reg",
        "blacklist"=list(blacklist)),
        nodes = bnlearn:::cache.structure(names(mb), arcs = arcs),
        arcs = arcs, blacklist=blacklist)
    constraints <- bnlearn:::arcs.to.be.added(res$arcs,
        nodes, whitelist = res$learning$blacklist)
    start <- empty.graph(nodes = nodes)

    hc_res <- hc(data, start=start, blacklist=constraints)
    return(hc_res)
}


.ncvregCV <- function(X, y, nFolds, pen) {
    pen <- strsplit(pen, "_") %>% vapply("[", 1, FUN.VALUE="a")
    fit <- cv.ncvreg(X, y, alpha=1, penalty=pen, nfolds=nFolds)
    included_var_index <- as.numeric(which(fit$fit$beta[, fit$min][-1]!=0))
    included_var_index
}

.glmnetBIC <- function(X, y, nFolds=NULL, pen=NULL) {
    ## This corresponds to L1MB
    fit <- glmnet(X, y, alpha=1, family="gaussian")
    tLL <- fit$nulldev - fit$nulldev * (1-fit$dev.ratio)
    k <- fit$df
    n <- nobs(fit)
    BIC <- log(n)*k - tLL
    included_var_index <- which(as.numeric(fit$beta[,which.min(BIC)])!=0)
    included_var_index
}

.plassoBIC <- function(X, y, nFolds=NULL, pen=NULL) {
    fit <- plasso_fit(X %>% as.matrix(), y, lambda=0.1, maxIter=100, gamma=0.5, eps=1e-6)
    ## The last coefficient is intercept
    included_var_index <- which(fit[1:(length(fit)-1)]!=0)
    included_var_index
}

.glmnetCV <- function(X, y, nFolds=5, pen=NULL) {
    fit <- cv.glmnet(X %>% as.matrix(), y, alpha=1, family="gaussian", nfolds=nFolds)
    numcoef <- coef(fit, s="lambda.min")[,1]
    numcoef <- numcoef[2:length(numcoef)]
    included_var_index <- which(numcoef!=0)
    included_var_index
}

.L0CV <- function(X, y, nFolds=5, pen=NULL) {
    fit=L0Learn::L0Learn.cvfit(x=X %>% as.matrix(), y=y, penalty="L0",
        nFolds=nFolds, algorithm="CD")
    lambdaIndex <- which.min(fit$cvMeans[[1]])
    cand_lambda <- fit$fit$lambda[[1]][lambdaIndex]
    L0coef <- coef(fit$fit, lambda=cand_lambda)
    vars <- as.numeric(L0coef)[2:length(L0coef)]
    included_var_index <- which(vars!=0)
    return(included_var_index)
}


.L0LXCV <- function(X, y, nFolds=5, pen=NULL) {
    pen <- strsplit(pen, "_") %>% vapply("[", 1, FUN.VALUE="a")
    fit=L0Learn::L0Learn.cvfit(x=X, y=y, penalty=pen, nFolds=nFolds, algorithm="CD")
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


