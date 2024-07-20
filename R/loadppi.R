# cache <- new.env()

#' loadppi
#' 
#' The igraph represenatation of the STRING database.
#' The file was downloaded from https://string-db.org/cgi/download, with the organism
#' corresponding to Homo sapiens or Mus musculus. The network was subset by the combined scores (above 900)
#' and attached to the package.
#' 
#' Szklarczyk D, Gable AL, Lyon D, Junge A, Wyder S, Huerta-Cepas J, Simonovic M, Doncheva NT, Morris JH, Bork P, Jensen LJ, Mering CV. STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res. 2019 Jan 8;47(D1):D607-D613. doi: 10.1093/nar/gky1131. PMID: 30476243; PMCID: PMC6323986.
#' 
#' @export
loadppi <- function(org="mm", database="string") {
    # if(is.null(cache$table)) {
    if (database=="string") {
      tb = readRDS(system.file("extdata", "ppi_mm_900.rds", package = "scstruc"))
    } else {
      tb = readRDS(system.file("extdata", "ppi_mm_biogrid.rds", package = "scstruc"))
    }
    # tb = cache$table
    tb
}


#' intersectPpi
#' 
#' Examine intersection of PPI using output of `markerCoefs`
#' 
#' @param edge_names edge names in `GeneA->GeneB` style
#' @param org organism name (mm, or hsa)
#' @export
intersectPpi <- function(edge_names, org="mm", database="string") {
  if (!is.matrix(edge_names)) {
    enn <- length(edge_names)
    importantGraph <- do.call(rbind,
                            edge_names %>% strsplit("->")) %>%
    graph_from_edgelist(directed = FALSE)    
  } else {
    enn <- dim(edge_names)[1]
    importantGraph <- graph_from_edgelist(edge_names, directed=FALSE)
  }

    ppis <- loadppi(org=org, database=database)
    ints <- igraph::intersection(ppis, importantGraph)

    inten <- as_edgelist(ints) %>% dim() %>% purrr::pluck(1)

    cat(inten / enn, "\n")
    return(ints)
}