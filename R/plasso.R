#' @useDynLib scstruc, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL
#' plasso_fit using RcppArmadillo is faster.
#' @noRd
plasso.fit <- function(X, y, lambdas=NULL, lambdas.length=10, P=NULL,
	eps=1e-6, maxIter=100, gamma=0.5, lr=1e-6, tol=1e-5, mu=1e-2) {
	## Generate lambdas using sparsebnUtils
    if (!is.matrix(X)) {stop("X must be matrix")}
	rlam <- 1e-2
	nn <- nrow(X)
    if (!requireNamespace("sparsebnUtils")) {
        stop("Needs installation of sparsebnUtils")
    } else {
        requireNamespace("sparsebnUtils")
    }
	lambdas <- sparsebnUtils::generate.lambdas(lambda.max = sqrt(nn),
       lambdas.ratio = rlam,
       lambdas.length = as.integer(lambdas.length),
       scale = "log")
	all_res <- lapply(lambdas, function(l) {
		fit <- plasso_fit(X, y,
			maxIter=maxIter, lambda=l, gamma=gamma, eps=eps)
		tmp_coef <- fit[1:(length(fit)-1)]
		pred <- plasso.predict(X, fit)
		## Remove intercept
		list("coef"=tmp_coef,
			"suppsize"=length(which(tmp_coef!=0)),
			"pred"=pred)
	})
	names(all_res) <- as.character(lambdas)
	return(all_res)
}


#' plasso.fit.single
#' 
#' @description reimplementation of PrecisionLasso.
#' Currently, only the continuous response is supported.
#' Automatic determination of gamma will be supported.
#' RcppArmadillo version is faster (plasso_fit).
#'  
#' @noRd
plasso.fit.single <- function(X, y, lambda=1, eps=1e-6, P=NULL,
    maxIter=100, gamma=0.5, lr=1e-6, tol=1e-5, mu=1e-2) {
    if (!is.matrix(X)) {stop("X must be matrix")}
    X0 <- cbind(rep(1, length(y)))
    X <- cbind(X, X0)
    P <- X
    n <- dim(X)[1]
    p <- dim(X)[2]
    k_max <- 0
    w <- rep(0, p)
    b <- t(X) %*% y
    G <- t(X) %*% X

    muList <- 10**(2 -8 * (seq_len(maxIter)-1) / maxIter)
    for (it in seq_len(maxIter)) {
        mu <- muList[it]
        mu <- lambda
        w[is.na(w)] <- 0
        w[abs(w) < mu] <- mu
        W <- P * t(replicate(dim(P)[1], w))
        eig <- eigen(W %*% t(W))
        s1 <- eig$values
        U1 <- eig$vectors
        s1[abs(s1) < mu] <- mu
        s1 <- abs(s1)

        W <- P * t(replicate(dim(P)[1], 1/w))
        eig2 <- eigen( W %*% t(W))
        s2 <- eig2$values
        U2 <- eig2$vectors
        s2[abs(s2) < mu] <- mu
        s2 <- abs(1/s2)

        s1 <- sqrt(s1)
        s2 <- sqrt(s2)
                
        U1 <- t(P) %*% U1
        D1 <- apply(U1 * U1 * t(replicate(dim(U1)[1], 1/s1)), 1, sum)
        maxi <- max(D1)
        mini <- min(D1)
        D1 <- (D1-mini) / (maxi-mini)

        U2 <- t(P) %*% U2
        D2 <- apply(U2 * U2 * t(replicate(dim(U2)[1], 1/s2)), 1, sum)
        maxi <- max(D2)
        mini <- min(D2)
        D2 <- (D2-mini) / (maxi-mini)

        D <- lambda * (gamma*D1 + (1-gamma)*D2)
        w_old <- w

        w <- solve(G + diag(D), b)
        w <- as.numeric(w[,1])
        if (sqrt(sum((w-w_old)**2)) < eps) {break}
    }
    w[is.na(w)] <- 0
    w[abs(w) < mu] <- 0
    return(w)
}

#' plasso.predict
#' prediction function for precision lasso
#' @param X prediction vector
#' @param w coefficient
#' @export
plasso.predict <- function(X, w) {
    X0 = rep(1, nrow(X))
    X = cbind(X, X0)
    X %*% w
}