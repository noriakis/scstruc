#' evaluate_metrics
#' 
#' s0: edge number in original graph
#' p: node number in original graph
#' N: sample number
#' 
#' @export
evaluate_metrics <- function(fitted, N, pens=c("glmnet_CV"),
         orig=TRUE, alphas=c(0.001, 0.005, 0.01, 0.05),
         ccdr=TRUE, return_net=FALSE) {
  rawnet <- as.bn(bn_fit_to_igraph(fitted))
  input <- rbn(fitted, N)
  alls <- lapply(pens, function(p) {
    cat(p,"")
    s <- Sys.time()
    net <- skeleton.reg(input, p, s="lambda.min")
    e <- Sys.time()
    tim <- as.numeric(e-s, unit="secs")
    cat(tim, "\n")
    list(net, tim)
  })
  names(alls) <- pens
  
  if (orig) {
    for (al in alphas) {
      s <- Sys.time()
      net <- mmhc(input, restrict.args=list("alpha"=al))
      e <- Sys.time()
      tim <- as.numeric(e-s, unit="secs")
      cat(al, tim, "\n")
      alls[[paste0("orig_",al)]] <- list(net, tim)
    }
  }
  if (ccdr) {
    s <- Sys.time()
    ccdr <- ccdr.run(sparsebnUtils::sparsebnData(input, type="continuous"),
                     lambdas.length=20)
    e <- Sys.time()
    for (i in ccdr) {
      alls[[paste0("ccdr_",round(i$lambda,3))]] <- list(sparsebnUtils::to_bn(i)$edges,
                                               as.numeric(e-s, unit="secs"))
    }  
  } 
  res <- do.call(rbind, lapply(names(alls), function(x) {
    cur_net <- alls[[x]][[1]]
    tim <- alls[[x]][[2]]
    s0 <- dim(rawnet$arcs)[1]
    comp <- bnlearn::compare(cur_net, rawnet)
    tp <- comp$tp; fp <- comp$fp; fn <- comp$fn
    pre <- tp/(tp+fp); rec <- tp/(tp+fn)
    fv <- 2*(pre*rec)/(pre+rec)
    c(
      x, s0,
      KL(bn.fit(cur_net, input), orig.net),
      BIC(cur_net, input),
      bnlearn::shd(cur_net, rawnet),
      tp, fp, fn, tp/s0, pre, rec, fv , tim
    )
  })) %>% data.frame()
  res <- res %>%  `colnames<-`(c("algo","s0","KL","BIC","SHD","TP","FP","FN","TPR","Precision","Recall","Fvalue","time")) %>%
    mutate_at(2:ncol(res), as.numeric)
  res <- res %>% mutate(BICnorm = (BIC - min(BIC)) / (max(BIC) - min(BIC)))
  res$N <- N
  res$p <- dim(input)[2]

  netList <- lapply(alls, function(x) x[[1]])
  names(netList) <- names(alls)
  if (return_net) {
     return(list(res, netList))    
  } else {
     return(res)
  }
}
