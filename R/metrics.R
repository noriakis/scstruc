#' @title metrics
#' 
#' @description Calculate and return metrics
#' 
#' @param reference reference bn object
#' @param inferred named list of inferred networks. Each element may be a bn object,
#' adjacency matrix, or edge list with score column.
#' @param sid_sym symmetric version of sid will be computed
#' @param SID.cran use implementation in package SID
#' @param data if provided, calculate KL, BIC, and AUPRC when possible
#' @return data.frame storing evaluation results
#' @export
#' @examples
#' library(bnlearn)
#' data(gaussian.test)
#' ref <- hc(head(gaussian.test))
#' metrics(ref, list("ref1"=ref,"ref2"=ref))
metrics <- function(reference, inferred, sid_sym=FALSE, SID.cran=FALSE, data=NULL) {
    if (is.null(names(inferred))) {
        stop("The inferred list should be named.")
    }
    infer_mode <- function(x) {
        if (is.matrix(x)) {
            if (nrow(x) == ncol(x) && isTRUE(all.equal(x, t(x), check.attributes = FALSE))) {
                return("undirected")
            }
            return("directed")
        }
        if (is.data.frame(x) && all(c("from", "to") %in% colnames(x))) {
            if (nrow(x) == 0) {
                return("directed")
            }
            key <- paste(x$from, x$to, sep = "\r")
            rev_key <- paste(x$to, x$from, sep = "\r")
            if (all(key %in% rev_key | rev_key %in% key)) {
                return("undirected")
            }
        }
        "directed"
    }
    count_edges <- function(x, mode) {
        if (inherits(x, "bn") || inherits(x, "bn.fit")) {
            amat <- bnlearn::amat(x)
            diag(amat) <- 0
            return(sum(amat != 0, na.rm = TRUE))
        }
        if (is.matrix(x)) {
            mat <- suppressWarnings(matrix(as.numeric(x),
                nrow = nrow(x), dimnames = dimnames(x)))
            diag(mat) <- 0
            if (mode == "undirected") {
                return(sum(mat != 0, na.rm = TRUE) / 2)
            }
            return(sum(mat != 0, na.rm = TRUE))
        }
        if (is.data.frame(x) && all(c("from", "to") %in% colnames(x))) {
            x <- x[x$from != x$to, , drop = FALSE]
            if (mode == "undirected") {
                keys <- paste(pmin(x$from, x$to), pmax(x$from, x$to), sep = "\r")
                return(length(unique(keys)))
            }
            return(nrow(unique(x[c("from", "to")])))
        }
        NA
    }
    res <- do.call(rbind, lapply(names(inferred), function(x) {
        cur_net <- inferred[[x]]
        s0 <- dim(reference$arcs)[1]
        rn <- length(reference$nodes)

        if (inherits(cur_net, "bn") || inherits(cur_net, "bn.fit")) {
            edges <- dim(cur_net$arcs)[1]
            comp <- bnlearn::compare(target=reference, current=cur_net)
            tp <- comp$tp; fp <- comp$fp; fn <- comp$fn
            pre <- tp/(tp+fp); rec <- tp/(tp+fn)
            fv <- 2*(pre*rec)/(pre+rec)
            auprc.val <- tryCatch({
                as.numeric(calc.auprc(reference, bnlearn::amat(cur_net))$auprc$.estimate[[1]])
            }, error = function(e) {
                NA
            })

            if (SID.cran) {
                sid.val <- SID.sid(reference, cur_net, sid_sym)
            } else {
                if (sid_sym) {
                    sid.val <- (bnlearn::sid(reference, cur_net) + bnlearn::sid(cur_net, reference))/2
                } else {
                    sid.val <- bnlearn::sid(cur_net, reference)
                }
            }
            inn <- length(cur_net$nodes)

            if (!is.null(data)) {
                kl <- bnlearn::KL(bnlearn::bn.fit(cur_net, data),
                	bnlearn::bn.fit(reference,data))
                bic <- BIC(cur_net, data)
            } else {
                kl <- NA
                bic <- NA
            }
            shd.val <- bnlearn::shd(cur_net, reference)
        } else {
            mode <- infer_mode(cur_net)
            edges <- count_edges(cur_net, mode)
            auprc.val <- tryCatch({
                as.numeric(calc.auprc(reference, cur_net, mode = mode)$auprc$.estimate[[1]])
            }, error = function(e) {
                NA
            })
            if (is.matrix(cur_net)) {
                inn <- ncol(cur_net)
            } else {
                inferred_nodes <- unique(c(cur_net$from, cur_net$to))
                inn <- length(inferred_nodes)
            }
            tp <- NA
            fp <- NA
            fn <- NA
            pre <- NA
            rec <- NA
            fv <- NA
            sid.val <- NA
            kl <- NA
            bic <- NA
            shd.val <- NA
        }

        c(
            x, rn, inn, s0, edges,
            shd.val,
            tp, fp, fn, tp/s0, pre, rec, fv, auprc.val, sid.val, kl, bic
        )
    })) %>% data.frame()
    res <- res %>%  `colnames<-`(c("algo","referenceNode", "InferenceNode",
        "s0","edges","SHD","TP",
        "FP","FN","TPR","Precision","Recall","F1","AUPRC","SID","KL","BIC")) %>%
        mutate_at(2:ncol(res), as.numeric)
    # res <- res %>% mutate(BICnorm = (BIC - min(BIC)) / (max(BIC) - min(BIC)))
    return(res)
}
