#' getKEGGPathwayGenes
#' 
#' 
#' @importFrom AnnotationDbi mapIds
#' @export
getKEGGPathwayGenes <- function(pid, org="mmu", orgDb=org.Mm.eg.db, type="SYMBOL") {
    path_eg  <- keggLink("pathway", org) %>% 
      tibble(pathway = ., eg = sub(paste0(org,":"), "", names(.))) %>%
      mutate(
        symbol = mapIds(orgDb, eg, type, "ENTREZID")
      )
    inc <- as.character(path_eg %>% filter(pathway == paste0("path:", pid)) %>% pull(symbol))
    inc
}


#' getGOGenes
#' 
#' get the vector of genes from GO term
#' 
#' @export
#' @param GO GO ID
#' @param orgDb AnnotationDbi that has "GOALL" slot
#' @param type return type, default to SYMBOL
getGOGenes <- function(GO, orgDb=org.Mm.eg.db, type="SYMBOL") {
  retrieved <- AnnotationDbi::select(org.Mm.eg.db,
     keytype="GOALL", keys=GO, columns="SYMBOL")
  retrieved$SYMBOL
}
