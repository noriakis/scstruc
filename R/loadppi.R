#' @title loadppi
#' @description The igraph represenatation of the STRING database.
#' The file was downloaded from https://string-db.org/cgi/download, with the organism
#' corresponding to Homo sapiens or Mus musculus. The network was subset by the combined scores (above 900)
#' and attached to the package.
#' 
#' Szklarczyk D, Gable AL, Lyon D, Junge A, Wyder S, Huerta-Cepas J, Simonovic M, Doncheva NT, Morris JH, Bork P, Jensen LJ, Mering CV. STRING v11: protein-protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. Nucleic Acids Res. 2019 Jan 8;47(D1):D607-D613. doi: 10.1093/nar/gky1131. PMID: 30476243; PMCID: PMC6323986.
#' 
#' @export
#' @param org organism, mm or hsa
#' @param database string or biogrid
#' @examples
#' db <- loadppi()
loadppi <- function(org="mm", database="string") {
    # if(is.null(cache$table)) {
    if (database=="string") {
      if (org == "mm") {
        tb = readRDS(system.file("extdata", "ppi_mm_900.rds", package = "scstruc"))
      } else if (org == "hsa") {
        tb = readRDS(system.file("extdata", "ppi_hs_900.rds", package = "scstruc"))
      } else {
        stop("mm or hsa should be specified to database")
      }
    } else if (database=="biogrid") {
      if (org == "mm") {
        tb = readRDS(system.file("extdata", "ppi_mm_biogrid.rds", package = "scstruc"))

      } else if (org=="hsa") {
        return(NULL)
      } else {
        stop("mm or hsa should be specified to database")
      }
    } else {
      stop("`database` argument should be string or biogrid")
    }
    tb
}

#' @title intersectPpi
#' 
#' @description
#' Examine intersection of PPI.
#' Named vector of edge number, intersection number and proportion will be returned.
#' 
#' @param edge_names edge names in `GeneA->GeneB` style character vector or igraph or bn
#' @param org organism name (mm, or hsa)
#' @param return_net returns the intersected igraph object
#' @export
#' @examples
#' intersectPpi(c("Gene1->Gene2"))
intersectPpi <- function(edge_names, org="mm", return_net=FALSE) {
  if (is.vector(edge_names)) {
      enn <- length(edge_names)
      importantGraph <- do.call(rbind,
                            edge_names %>% unique() %>% strsplit("->")) %>%
      igraph::graph_from_edgelist(directed = FALSE)
  } else if ("bn" %in% class(edge_names)) {
    ## If BN, considering the undirected edges, simplify
    g <- bnlearn::as.igraph(edge_names) %>% 
        igraph::as.undirected() %>% igraph::simplify()
    enn <- length(igraph::E(g))
    importantGraph <- g
  } else if (igraph::is.igraph(edge_names)) {
    ## If igraph, simplify
    importantGraph <- igraph::as.undirected(edge_names) %>% igraph::simplify()
    enn <- igraph::E(importantGraph) %>% length()
  } else if (is.matrix(edge_names)) {
    ## User-specification warning
    warning("Matrix will be taken care of as is, so multiple rows with same relationships (undirected edges) will be counted as two.")
    enn <- dim(edge_names)[1]
    importantGraph <- igraph::graph_from_edgelist(edge_names, directed=FALSE) %>% igraph::simplify()
  }

    ppis <- loadppi(org=org, database="string")
    ints <- igraph::intersection(ppis, importantGraph)

    ppis.biogrid <- loadppi(org=org, database="biogrid")
    ints.biogrid <- igraph::intersection(ppis.biogrid, importantGraph)

    if (return_net) {
      return(list("string"=ints, "biogrid"=ints.biogrid))
    }

    inten <- dim(igraph::as_edgelist(ints))[1]
    inten.biogrid <- dim(igraph::as_edgelist(ints.biogrid))[1]

    # cat(inten / enn, "\n")
    vec <- c(enn, inten, inten / enn, "string")
    names(vec) <- c("edge_number", "intersect", "proportion", "database")

    vec2 <- c(enn, inten.biogrid, inten.biogrid / enn, "biogrid")
    names(vec2) <- c("edge_number", "intersect", "proportion", "database")

    df <- data.frame(rbind(vec, vec2))
    df$edge_number <- as.numeric(df$edge_number)
    df$intersect <- as.numeric(df$intersect)
    df$proportion <- as.numeric(df$proportion)
    row.names(df) <- seq_len(nrow(df))
    return(df)
}

