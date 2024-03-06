
#' @import AnnotatiobDbi mapIds
returnKEGGPathSymbol <- function(pid, org="mmu", orgDb=org.Mm.eg.db) {
    path_eg  <- keggLink("pathway", org) %>% 
      tibble(pathway = ., eg = sub(paste0(org,":"), "", names(.))) %>%
      mutate(
        symbol = mapIds(orgDb, eg, "SYMBOL", "ENTREZID")
      )

    inc <- as.character(path_eg %>% filter(pathway == "path:mmu04150") %>% pull(symbol))
    inc
}

