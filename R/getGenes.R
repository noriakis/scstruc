#' @title getKEGGPathwayGenes
#' @description get genes related to KEGG PATHWAY
#' @param pid pathway id in KEGG
#' @param org organism name
#' @param orgDb organism database
#' @param type returned ID type of genes
#' @importFrom tibble tibble
#' @export
getKEGGPathwayGenes <- function(pid, org="mmu", orgDb=NULL, type="SYMBOL") {
    if (!requireNamespace("AnnotationDbi")) {
      stop("Needs AnnotationDbi installation")
    }
    if (is.null(orgDb)) {
        stop("Please specify annotation database e.g. org.Mm.eg.db")
    }
    path_eg  <- KEGGREST::keggLink("pathway", org)
    path_eg <-  tibble::tibble(pathway = path_eg,
        eg = sub(paste0(org,":"), "", names(path_eg))) %>%
      mutate(
        symbol = AnnotationDbi::mapIds(orgDb, .data$eg, type, "ENTREZID")
      )
    inc <- as.character(path_eg %>%
        dplyr::filter(.data$pathway == paste0("path:", pid)) %>% pull(.data$symbol))
    inc
}

#' @title getGOGenes
#' @description get the vector of genes related to GO term
#' @export
#' @param GO GO ID
#' @param orgDb AnnotationDbi that has "GOALL" slot
#' @param type return type, default to SYMBOL
getGOGenes <- function(GO, orgDb=NULL, type="SYMBOL") {
  if (!requireNamespace("AnnotationDbi")) {
    stop("Needs AnnotationDbi installation")
  }
    if (is.null(orgDb)) {
        stop("Please specify annotation database e.g. org.Mm.eg.db")
    }
    retrieved <- AnnotationDbi::select(orgDb,
        keytype="GOALL", keys=GO, columns="SYMBOL")
    retrieved$SYMBOL
}
