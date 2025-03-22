#' @title metricsFromFitted
#' 
#' @description 
#' Evaluate various metrics comparing reference network and inferred network.
#' The inference will be performed based on randomly generated data from 
#' fitted Bayesian network (`bn.fit` object). For the already calculated networks,
#' please use `metrics()` function.
#' 
#' The description of metrics is as follows:
#' 
#' s0: edge number in original graph
#' edges: edge number in inferred graph
#' p: node number in original graph
#' N: sample number
#' SHD: structural hamming distance
#' SID: structural intervention distance (if specified)
#' time: time needed in structure learning (sec)
#' 
#' Time will be calculated by differences in Sys.time() and will be displayed after
#' algorithms' name. As an external software, GENIE3 will be tested if specified.
#' 
#' @param fitted reference network that is parameter fitted
#' @param N number of samples to be drawn from `fitted` by `rbn` function.
#' @param mmhc perform comparison with MMHC implemented in bnlearn with
#' alphas specified in `alphas` argument.
#' @param alphas alpha to be used in the restricting phase of MMHC.
#' @param algos algorithm names
#' (For CCDr, hurdle, and MMHC, please specify other arguments)
#' @param hurdle perform inference based on Hurdle model, default to FALSE.
#' Some scoring functions used in the function cannot calculate scores
#' if the feature abundances are all negative. In this case, these scores will be omitted.
#' @param ccdr perform comparison with CCDr algorithm with lambda length specified in 
#' `lambdas.length` argument.
#' @param ppi if specified, include metrics of comparison with PPI.
#' The network node name should be symbol for calculating the intersection.
#' Also, `database` and `org` argument should be specified based on the node name.
#' @param org passed to `intersectPpi` function
#' @param genie test genie3 for the comparison
#' @param genie.threshold genie3 threshold value. The weight has no statistical meanings,
#' so the value should be chosen carefully.
#' @param sid compute SID
#' @param sid_sym compute symmetrized version of SID
#' @param return_net return list of whole BN
#' @param return_data return data
#' @param SID.cran use package SID for calculation of SID
#' @param ges perform ges evaluation
#' @param lingam perform lingam evaluation
#' @param lambdas.length lambda length in CCDr algorithm
#' @param algorithm.args arguments to be passed to algorithms,
#' must be the same length as `algos`.
#' @export
metricsFromFitted <- function(fitted, N, algos=c("glmnet_CV"),
    mmhc=TRUE, alphas=c(0.001, 0.005, 0.01, 0.05), ppi=FALSE,
    org="mm", return_data=FALSE, SID.cran=FALSE,
    algorithm.args=list(), hurdle=FALSE, sid_sym=FALSE,
    ccdr=TRUE, return_net=FALSE, lingam=FALSE,
    genie=FALSE, genie.threshold=seq(0,1,0.1),
    lambdas.length=10, sid=FALSE, ges=FALSE) {

    if ("ccdr" %in% algos) {
      stop("Cannot specify CCDr in this argument, please choose `ccdr` argument as TRUE")
    }
    if ("ccdr.boot" %in% algos) {
      stop("Cannot specify CCDr in this argument, please choose `ccdr` argument as TRUE")
    }
    if ("Hurdle" %in% algos) {
      stop("Cannot specify Hurdle in this argument, please choose `hurdle` argument as TRUE")
    }

    # if (length(algorithm.args)!=0) {
    #     if (length(algorithm.args)!=length(algos)){
    #         stop("Length of algos and algorithm.args does not match")
    #     }
    # }

    rawnet <- bnlearn::as.bn(bn_fit_to_igraph(fitted))
    input <- bnlearn::rbn(fitted, N)
    alls <- lapply(algos, function(p) {
        cat_subtle(p," ")
        s <- Sys.time()
        net <- skeleton.reg(input, p, s="lambda.min")
        e <- Sys.time()
        tim <- as.numeric(e-s, unit="secs")
        cat_subtle(tim, "\n")
        list(net, tim)
    })
    names(alls) <- algos

    if (hurdle) {
        ###
        ## Two types of scores will be tested
        ## Hurdle object can be passed to .Hurdle, but
        ## for evaluation of inferencetime,
        ## raw calculation will be performed.
        ###

        s <- Sys.time()
        net <- .Hurdle(input, score=NULL)
        e <- Sys.time()
        tim <- as.numeric(e-s, unit="secs")
        cat_subtle("Hurdle BIC ", tim, "\n")
        alls[[paste0("Hurdle_BIC")]] <- list(net, tim)

        s <- Sys.time()
        net <- .Hurdle(input, score=hurdle.bic)
        e <- Sys.time()
        tim <- as.numeric(e-s, unit="secs")
        cat_subtle("Hurdle zBIC ", tim, "\n")
        alls[[paste0("Hurdle_zBIC")]] <- list(net, tim)

    }

    if (genie) {
    	if (!requireNamespace("GENIE3")) {
    		stop("Please install GENIE3")
    	}
        genie.input <- input %>% t()
        s <- Sys.time()
        gra <- GENIE3::GENIE3(genie.input)
        e <- Sys.time()
        nets <- lapply(genie.threshold, function(th) {
            tmp <- gra
            tmp[tmp < th] <- 0
            tmp.net <- igraph::graph_from_adjacency_matrix(tmp, weighted = TRUE)
            if (igraph::is.dag(tmp.net)) {
                return(bnlearn::as.bn(tmp.net))
            } else {
                return(NA)
            }
        })
        names(nets) <- genie.threshold
        if (sum(unlist(lapply(nets, function(net) {is.na(net)}))) == length(genie.threshold)) {
            cat_subtle("GENIE3 All threshold results in non-DAG\n")
        } else {
            tim <- as.numeric(e-s, unit="secs")
            cat_subtle("GENIE3 ", tim, "\n")
            for (nn in names(nets)) {
                genie.net <- nets[[nn]]
                if (is(genie.net, "bn")) {
                    alls[[paste0("GENIE3_",nn)]] <- list(genie.net, tim)
                }
            }
        }
    }

    if (ges) {
        tryCatch({
            s <- Sys.time()
            net <- ges.pcalg(input)
            e <- Sys.time()
            tim <- as.numeric(e-s, unit="secs")
            cat_subtle("GES ", tim, "\n")
            alls[[paste0("GES")]] <- list(net, tim)
            
        }, error= function(e) message(e))
    }
    
    if (lingam) {
        tryCatch({
            s <- Sys.time()
            net <- lingam.pcalg(input)
            e <- Sys.time()
            tim <- as.numeric(e-s, unit="secs")
            cat_subtle("LiNGAM ", tim, "\n")
            alls[[paste0("LiNGAM")]] <- list(net, tim)
        }, error = function(e) message(e))
    }


    if (mmhc) {
        for (al in alphas) {
            s <- Sys.time()
            net <- mmhc(input, restrict.args=list("alpha"=al))
            e <- Sys.time()
            tim <- as.numeric(e-s, unit="secs")
            cat_subtle("MMHC ", al, " ", tim, "\n")
            alls[[paste0("mmhc_",al)]] <- list(net, tim)
        }
    }
    if (ccdr) {
	    if (!requireNamespace("ccdrAlgorithm")) {
	        stop("Needs installation of ccdrAlgorithm")
	    } else {
	        requireNamespace("ccdrAlgorithm")
	    }
	    if (!requireNamespace("sparsebnUtils")) {
	        stop("Needs installation of sparsebnUtils")
	    } else {
	        requireNamespace("sparsebnUtils")
	    }
        s <- Sys.time()
        ccdr <- ccdrAlgorithm::ccdr.run(sparsebnUtils::sparsebnData(input, type="continuous"),
            lambdas.length=lambdas.length)
        e <- Sys.time()
        for (i in ccdr) {
            alls[[paste0("ccdr_",round(i$lambda,3))]] <- list(sparsebnUtils::to_bn(i)$edges,
                                               as.numeric(e-s, unit="secs"))
        }    
    } 
    cat_subtle("Network computing finished", "\n")

    res <- do.call(rbind, lapply(names(alls), function(x) {
        cur_net <- alls[[x]][[1]]
        tim <- alls[[x]][[2]]
        s0 <- dim(rawnet$arcs)[1]
        edges <- dim(cur_net$arcs)[1]
        comp <- bnlearn::compare(cur_net, rawnet)
        tp <- comp$tp; fp <- comp$fp; fn <- comp$fn
        pre <- tp/(tp+fp); rec <- tp/(tp+fn)
        fv <- 2*(pre*rec)/(pre+rec)

        if (ppi) {
            numppi <- intersectPpi(cur_net %>% as.igraph() %>% igraph::as_edgelist(),
                org=org)
        } else {
            numppi <- NA
        }
        if (sid) {
            if (SID.cran) {
                sid.val.sym <- SID.sid(rawnet, cur_net, sid_sym)
            } else {
                if (sid_sym) {
                    sid.val <- bnlearn::sid(cur_net, rawnet)
                    sid.val.2 <- bnlearn::sid(rawnet, cur_net)
                    sid.val.sym <- (sid.val + sid.val.2) / 2
                } else {
                    sid.val.sym <- bnlearn::sid(cur_net, rawnet)
                }                
            }
        } else {
            sid.val.sym <- NA
        }
        c(
            x, s0, edges,
            bnlearn::KL(bnlearn::bn.fit(cur_net, input), fitted),
            BIC(cur_net, input),
            bnlearn::shd(cur_net, rawnet),
            tp, fp, fn, tp/s0, pre, rec, fv , sid.val.sym, numppi, tim
        )
    })) %>% data.frame()
    res <- res %>%  `colnames<-`(c("algo","s0","edges","KL","BIC","SHD","TP",
        "FP","FN","TPR","Precision","Recall","F1","SID","PPI","time")) %>%
        mutate_at(2:ncol(res), as.numeric)
    res <- res %>% mutate(BICnorm = (BIC - min(BIC)) / (max(BIC) - min(BIC)))
    res$N <- N
    res$p <- dim(input)[2]

    netList <- lapply(alls, function(x) x[[1]])
    names(netList) <- names(alls)

    returns <- list()
    returns[["metrics"]] <- res
    if (return_net) {
        returns[["net"]] <- netList
    }
    if (return_data) {
        returns[["data"]] <- input
    }
    return(returns)
}
