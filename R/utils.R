#' @noRd
load.grnboost2 <- function(filename, minmax=TRUE,
                           thresholds=seq(0, 1, 0.1)) {
    minmaxFn <- function(x, na.rm = TRUE) {
        return((x- min(x)) /(max(x)-min(x)))
    }
    df <- read.table(filename, sep="\t", header=1)
    if (minmax) {
        df$importance <- minmaxFn(df$importance)
    } else {
    }
    g <- igraph::graph_from_data_frame(df[, c("TF","target","importance")])
    adj <- as_adj(g, attr="importance")
    lapply(thresholds, function(th) {
        adj[adj<th] <- 0
        tmpg <- igraph::graph_from_adjacency_matrix(adj, weighted=TRUE)
        
        if (igraph::is_dag(tmpg)) {
            if (length(V(tmpg))!=0) {
                return(bnlearn::as.bn(tmpg))
            } else {
                return(NA)
            }
        } else {
            return(NA)
        }
    }) %>% setNames(paste0("GRNBoost2_",thresholds))
}

#' @noRd
bn.recovery.2 <- function(mb) {
    newmb <- list()
    nodes <- names(mb)
    for (node in nodes) {
        new <- list("mb"=c(), "nbr"=c())
        new[["mb"]] <- mb[[node]]$mb
        for (node2 in nodes) {
            c1 <- node2 %in% mb[[node]]$nbr
            c2 <- node %in% mb[[node2]]$nbr
            c3 <- c1 + c2
            if (c3 == 0) {} else if (c3 == 2) {new[["nbr"]] <- c(new[["nbr"]], node2)} else {}
        }
        newmb[[node]] <- new
    }
    return(newmb)
}


#' @noRd
arcs.to.be.added.2 <- function(arcs, nodes) {
    p <- length(nodes)
    mat <- matrix(1, nrow=p, ncol=p)
    row.names(mat) <- nodes
    colnames(mat) <- nodes
    for (nr in seq_len(nrow(arcs))) {
        tmp <- arcs[nr, ]
        a <- tmp[1]; b <- tmp[2]
        mat[a, b] <- 0
        mat[b, a] <- 0
    }
    igraph::as_edgelist(igraph::graph_from_adjacency_matrix(mat,diag = FALSE))
}


#' @noRd
bnlearn.sid.sym <- function(net1, net2) {
    sid1 <- bnlearn::sid(net1, net2)
    sid2 <- bnlearn::sid(net2, net1)
    (sid1 + sid2) / 2
}


#' @noRd
SID.sid <- function(trueBn, estBn, sym=FALSE) {
    rawadj.t <- as_adj(as.igraph(trueBn))
    rawadj.e <- as_adj(as.igraph(estBn))
    sid <- structIntervDist(rawadj.t, rawadj.e)$sid
    if (sym) {
        sid2 <- structIntervDist(rawadj.e, rawadj.t)$sid
        return((sid + sid2)/2)
    }
    return(sid)
}


#' @title bn.fit.hurdle
#' @description make a pseudo bn.fit object using continuous part of Hurdle model.
#' @details You should not use `rbn` or relevant functions for logic sampling 
#' to this object as  the logistic regression part of the model will 
#' not be taken into account.
#' @param x bn object
#' @param data data to be fitted
#' @param cdrAdjustment adjust for cellular detection rate based on `data`
#' @export
#' @examples
#' data(gaussian.test)
#' test.bn <- hc(gaussian.test)
#' bn.fit.hurdle(test.bn, gaussian.test)
bn.fit.hurdle <- function(x, data, cdrAdjustment=FALSE) {
    
    cdr <- rowSums(data>0)
    data[["cdr"]] <- scale(cdr)

    fits <- lapply(names(x$nodes), function(nn) {
        parents <- x$nodes[[nn]]$parents
        if (length(parents)==0) {
            if (cdrAdjustment) {
                mod <- paste0(nn, "~ cdr")
            } else {
                mod <- paste0(nn, " ~ 1")
            }
        } else {
            if (cdrAdjustment) {
                mod <- paste0(nn, "~", paste0(parents, collapse=" + "), " + cdr")
            } else {
                mod <- paste0(nn, "~", paste0(parents, collapse=" + "))
            }
        }
        fit <- MAST::zlm(as.formula(mod), data=data)
        res <- list()
        res[["node"]] <- nn
        res[["parents"]] <- parents
        res[["children"]] <- x$nodes[[nn]]$children
        coefs <- coef(fit$cont)
        coefs <- coefs[names(coefs)!="cdr"]
        res[["coefficients"]] <- coefs
        res[["residuals"]] <- fit$cont$residuals
        res[["sd"]] <- sd(fit$cont$residuals)
        res[["fitted.values"]] <- fit$cont$fitted.values
        ## Additionally storing fitted object
        res[["fitted"]] <- fit
        class(res) <- c("bn.fit.gnode")
        res
    }) %>% setNames(names(x$nodes))
    class(fits) <- c(class(fits), "bn.fit", "bn.fit.gnet")
    return(fits)
}



#' sel.genes.kegg.ora
#' select genes based on KEGG over-representation analysis
#' @param genes query genes
#' @param orgDb organism database object
#' @param keyType key type of genes
#' @param sel select the number of pathways, sorted by p-value
#' @param organism organism ID in KEGG, default to hsa
#' @export
sel.genes.kegg.ora <- function(genes, orgDb=org.Hs.eg.db, keyType="ENSEMBL", sel=NULL, organism="hsa") {
  if (!requireNamespace("AnnotationDbi")) {
    stop("Needs AnnotationDbi installation")
  }
  if (!requireNamespace("clusterProfiler")) {
    stop("Needs clusterProfiler installation")
  }
  input <- AnnotationDbi::mapIds(orgDb, genes, column="ENTREZID", keytype=keyType, multiVals="first") %>%
    as.character()
  input <- input[!is.na(input)]
  ora <- enrichKEGG(input, organism=organism)
  if (!is.null(sel)) {
    cat(ora@result[sel, ]$Description, "\n")
    ora <- ora@result[sel, ]$geneID %>% strsplit("/") %>% unlist()
    ora <- AnnotationDbi::mapIds(orgDb, ora, column=keyType, keytype="ENTREZID",
                                 multiVals="first") %>% as.character()
  }
  return(ora)
}

#' @noRd
cat_subtle <- function(...) cat(pillar::style_subtle(paste0(...)))

#' shdmat
#' @param netlist list of bn
#' @export
#' @examples
#' data(gaussian.test)
#' test <- head(gaussian.test, 10)
#' shdmat(list(hc(test), mmhc(test)))
shdmat <- function(netlist) {
    sapply(netlist, function(x) sapply(netlist, function(y) bnlearn::shd(x,y)))
}

#' sidmat
#' @param netlist list of bn
#' @param SID.cran use implementation in package SID
#' @export
#' @examples
#' data(gaussian.test)
#' test <- head(gaussian.test, 10)
#' sidmat(list(hc(test), mmhc(test)))
sidmat <- function(netlist, SID.cran=FALSE) {
   sapply(netlist, function(x) {
       sapply(netlist, function(y) {
            if (SID.cran) {
                val <- SID.sid(x, y)
            } else {
                val <- bnlearn::sid(x, y)
            }
           # rawadj.x <- as_adj(as.igraph(x))
           # rawadj.y <- as_adj(as.igraph(y))
           # sid <- SID::structIntervDist(rawadj.x, rawadj.y)$sid
           return(val)
       })
    })
}


#' @noRd
removeAllNegative <- function(input) {
    cat_subtle("Removing all negative genes if available\n")
    negg <- lapply(colnames(input), function(x) {
        tmp <- input[,x]
        all(tmp[tmp!=0] < 0)
    }) %>% unlist()
    if (any(negg)) {
        cat_subtle("  Found some genes\n")
        cat_subtle("  ", paste(colnames(input)[negg], collapse=", "), "\n")
        input <- input[, !negg]
    } else {
        cat_subtle("  None available\n")
    }
    return(input)
}

#' convert_to_tbl_graph
#' convert bn object to tbl_graph object
#' @param bn bn object
#' @return tbl_graph
#' @export
#' @examples
#' data(gaussian.test)
#' test <- head(gaussian.test)
#' convert_to_tbl_graph(hc(test))
convert_to_tbl_graph <- function(bn) {
    bn |> 
    bnlearn::as.igraph() |>
    tidygraph::as_tbl_graph()
}

#' 
#' bn_fit_to_igraph
#' 
#' This function coverts the bn.fit object to igraph object.
#' The nodes with no parents and children will not be included in the resulting igraph.
#' The coefficient value is stored in `coef` attribute of resulting graph.
#' 
#' @param graph bn.fit object
#' @param preserve preserve the node with no parents and children
#' @export
#' @examples
#' data(gaussian.test)
#' test <- head(gaussian.test)
#' struc <- hc(test)
#' bn_fit_to_igraph(bn.fit(struc, test))
bn_fit_to_igraph <- function(graph, preserve=FALSE) {
	tmp <- igraph::graph_from_data_frame(do.call(rbind, lapply(names(graph), function(i) {
	    tmp <- graph[[i]]
	    if (length(tmp$coefficients)>1) {
	        tmp <- data.frame(tmp$coefficients[2:length(tmp$coefficients)]) %>%
	        `colnames<-`(c("coef"))
	        tmp[["from"]] <- row.names(tmp)
	        tmp[["to"]] <- rep(i, nrow(tmp))
	        return(tmp[,c("from","to","coef")])
	    }
	})))
    if (!preserve) {
        return(tmp)
    }
    tmp <- tmp %>% as_tbl_graph()
    inName <- tmp %N>% pull(name)
    dinn <-setdiff(names(graph), inName)
    
    if (length(dinn)==0) {return(tmp %>% as.igraph())}

    tmp %>% bind_nodes(data.frame(name=dinn)) %>% as.igraph()
}

#' 
#' bn_fit_to_df
#' 
#' This function coverts the bn.fit object to data.frame
#' 
#' @param fitted bn.fit object
#' @export
#' @examples
#' data(gaussian.test)
#' test <- head(gaussian.test)
#' struc <- hc(test)
#' bn_fit_to_df(bn.fit(struc, test))
bn_fit_to_df <- function(fitted) {
    do.call(rbind, lapply(names(fitted), function(x) {
        tmp.n <- names(fitted[[x]]$coefficients)
        ps <- tmp.n[tmp.n!="(Intercept)"]
        do.call(rbind, lapply(ps, function(p) {
            c(p, x, fitted[[x]]$coefficients[p]) %>% setNames(c("from","to","coefficient"))
        })) %>% data.frame()
    }))
}

#' getCyclinGenes
#' @noRd
getCyclinGenes <- function(sce, id="Symbol", title=FALSE) {
	if (title) {
		cyclin.genes <- grep("^Ccn[abde]+", rowData(sce)[["Symbol"]])
	} else {
	    cyclin.genes <- grep("^CCN[ABDE]+", rowData(sce)[["Symbol"]])
	}

    cyclin.genes <- rownames(sce)[cyclin.genes]
    return(cyclin.genes)
}
