#' @title returnKEGGedges
#' @description return directed igraph object
#' @param pathID pathway ID
#' @param args arguments of ggkegg::pathway
#' @param bn if TRUE, try to convert to bn object, but if failed, return igraph
#' @param largestComponents obtain largest components in the graph
#' @param removeCycle remove cycle using feedback_arc_set()
#' @details The function uses ggkegg to obtain and parse KEGG PATHWAY information
#' using KEGG RESTful API.  There should be type `gene` in the node data of pathway
#' this returns the directed relationship described in the pathway.
#' The function tries to extract DAG from KEGG PATHWAY, and not necessarily succeed.
#' The first part of graphics name in the graph is taken as the representative gene name.
#' @export
#' @importFrom dplyr filter select
#' @importFrom igraph feedback_arc_set
getKEGGedges <- function(pathID, args=list(), largestComponents=TRUE, bn=TRUE,
  removeCycle=FALSE) {
    if (!requireNamespace("ggkegg")) {
        stop("Needs ggkegg. Please install from Bioconductor.")
    }
    args[["pid"]] <- pathID
    pway <- do.call(ggkegg::pathway, args)
    gs <- pway %N>% dplyr::filter(type=="gene")
    if (largestComponents) {
        comps <- gs %>% to_components()
        ind <- which.max(lapply(comps, function(x) {dim(x %N>% data.frame())[1]}))
        gs <- comps[[ind]]      
    }
    nname <- gs %N>% data.frame() %>% dplyr::pull(graphics_name) %>% strsplit(",") %>% sapply("[",1) %>% strsplit("\\.\\.\\.") %>% sapply("[", 1)
    ig <- gs %E>% mutate(fromn = nname[from], ton=nname[to]) %>% data.frame() %>% mutate(from=fromn, to=ton) %>%
        dplyr::select(-fromn) %>% dplyr::select(-ton) %>% igraph::graph_from_data_frame(directed=TRUE) %>% simplify()


    if (removeCycle) {
      ## Use feedback_arc_set to find minimum feedback
      es <- feedback_arc_set(ig)
      # es.dat <- do.call(rbind, as_ids(es) %>% strsplit("\\|"))
      if (length(es)==0) {

      } else {
        cat_subtle("Removing ", paste0(as_ids(es), collapse=", "), "\n")
        ig <- delete_edges(ig, es)
      }

    }
    if (bn) {
      as.bn(ig)
    } else {
      ig
    }

}
