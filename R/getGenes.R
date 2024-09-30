#' @title getKEGGPathwayGenes
#' @description get genes related to KEGG PATHWAY
#' @param pid pathway id in KEGG
#' @param org organism name
#' @param orgDb organism database
#' @param type returned ID type of genes
#' @export
getKEGGPathwayGenes <- function(pid, org="mmu", orgDb=org.Mm.eg.db, type="SYMBOL") {
    if (!requireNamespace("AnnotationDbi")) {
      stop("Needs AnnotationDbi installation")
    }
    path_eg  <- keggLink("pathway", org) %>% 
      tibble(pathway = ., eg = sub(paste0(org,":"), "", names(.))) %>%
      mutate(
        symbol = mapIds(orgDb, eg, type, "ENTREZID")
      )
    inc <- as.character(path_eg %>% dplyr::filter(.data$pathway == paste0("path:", pid)) %>% pull(symbol))
    inc
}

#' @title getGOGenes
#' @description get the vector of genes related to GO term
#' @export
#' @param GO GO ID
#' @param orgDb AnnotationDbi that has "GOALL" slot
#' @param type return type, default to SYMBOL
getGOGenes <- function(GO, orgDb=org.Mm.eg.db, type="SYMBOL") {
  if (!requireNamespace("AnnotationDbi")) {
    stop("Needs AnnotationDbi installation")
  }
  retrieved <- AnnotationDbi::select(org.Mm.eg.db,
     keytype="GOALL", keys=GO, columns="SYMBOL")
  retrieved$SYMBOL
}
