
#' sel.genes.kegg.ora
#' select genes based on KEGG over-representation analysis
#' @export
sel.genes.kegg.ora <- function(genes, orgDb=org.Hs.eg.db, keyType="ENSEMBL", sel=NULL, organism="hsa") {
  input <- AnnotationDbi::mapIds(orgDb, genes, column="ENTREZID", keytype=keyType, multiVals="first") %>%
    as.character()
  input <- input[!is.na(input)]
  ora <- clusterProfiler::enrichKEGG(input, organism=organism)
  if (!is.null(sel)) {
    cat(ora@result[sel, ]$Description, "\n")
    ora <- ora@result[sel, ]$geneID %>% strsplit("/") %>% unlist()
    ora <- AnnotationDbi::mapIds(orgDb, ora, column=keyType, keytype="ENTREZID",
                                 multiVals="first") %>% as.character()
  }
  return(ora)
}


#' shdmat
#' 
#' @export
shdmat <- function(netlist) {
    sapply(netlist, function(x) sapply(netlist, function(y) bnlearn::shd(x,y)))
}


#' convert_to_tbl_graph
#' convert bn object to tbl_graph object
#' @param bn bn object
#' @return tbl_graph
#' @export
convert_to_tbl_graph <- function(bn) {
    bn |> 
    bnlearn::as.igraph() |>
    tidygraph::as_tbl_graph()
}

#' 
#' bn_fit_to_igraph
#' 
#' This function coverts the bn.fit object to igraph object.
#' The nodes with no parents or nodes will not be included in the resulting igraph.
#' The coefficient value is stored in `coef` attribute of resulting graph.
#' 
#' @param graph bn.fit object
#' 
#' @export
bn_fit_to_igraph <- function(graph) {
	igraph::graph_from_data_frame(do.call(rbind, lapply(names(graph), function(i) {
	    tmp <- graph[[i]]
	    if (length(tmp$coefficients)>1) {
	        tmp <- data.frame(tmp$coefficients[2:length(tmp$coefficients)]) %>%
	        `colnames<-`(c("coef"))
	        tmp[["from"]] <- row.names(tmp)
	        tmp[["to"]] <- rep(i, nrow(tmp))
	        return(tmp[,c("from","to","coef")])
	    }
	})))
}


#' getCyclinGenes
#' @export
getCyclinGenes <- function(sce, id="Symbol", title=FALSE) {
	if (title) {
		cyclin.genes <- grep("^Ccn[abde]+", rowData(sce)[["Symbol"]])
	} else {
	    cyclin.genes <- grep("^CCN[ABDE]+", rowData(sce)[["Symbol"]])
	}

    cyclin.genes <- rownames(sce)[cyclin.genes]
    return(cyclin.genes)
}
