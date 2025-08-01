#' @title metrics
#' 
#' @description Calculate and return metrics
#' 
#' @param reference reference bn object
#' @param inferred named list of bn object
#' @param sid_sym symmetric version of sid will be computed
#' @param SID.cran use implementation in package SID
#' @param data if provided, calculate KL
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
    res <- do.call(rbind, lapply(names(inferred), function(x) {
        cur_net <- inferred[[x]]
        s0 <- dim(reference$arcs)[1]
        edges <- dim(cur_net$arcs)[1]
        comp <- bnlearn::compare(target=reference, current=cur_net)

        tp <- comp$tp; fp <- comp$fp; fn <- comp$fn
        pre <- tp/(tp+fp); rec <- tp/(tp+fn)
        fv <- 2*(pre*rec)/(pre+rec)

        if (SID.cran) {
            sid.val <- SID.sid(reference, cur_net, sid_sym)
        } else {
            if (sid_sym) {
                sid.val <- (bnlearn::sid(reference, cur_net) + bnlearn::sid(cur_net, reference))/2
            } else {
                sid.val <- bnlearn::sid(cur_net, reference)
            }
        }
        rn <- length(reference$nodes)
        inn <- length(cur_net$nodes)

        if (!is.null(data)) {
            kl <- bnlearn::KL(bnlearn::bn.fit(cur_net, data),
            	bnlearn::bn.fit(reference,data))
            bic <- BIC(cur_net, data)
        } else {
            kl <- NA
            bic <- NA
        }

        c(
            x, rn, inn, s0, edges,
            # KL(bn.fit(cur_net, input), fitted),
            # BIC(cur_net, input),
            bnlearn::shd(cur_net, reference),
            tp, fp, fn, tp/s0, pre, rec, fv, sid.val, kl, bic
        )
    })) %>% data.frame()
    res <- res %>%  `colnames<-`(c("algo","referenceNode", "InferenceNode",
        "s0","edges","SHD","TP",
        "FP","FN","TPR","Precision","Recall","F1","SID","KL","BIC")) %>%
        mutate_at(2:ncol(res), as.numeric)
    # res <- res %>% mutate(BICnorm = (BIC - min(BIC)) / (max(BIC) - min(BIC)))
    return(res)
}