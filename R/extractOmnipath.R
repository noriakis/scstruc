#' @title extractDAGfromOmnipath
#' This function takes the pathway name in SignaLink and
#' 1. Subset the pathway-related genes in dorothea interaction and remove duplicates
#' 2. Identify and remove minimal feedback arc set (make DAG)
#' 3. Takes the largest components in the graph
#' @param genes gene symbols used to subset the network
#' @param mfas force DAG
#' @param pathway_name pathway name in SignaLink
#' @param target_source whether to subset by source/target pathway or genes using
#' `AND` or `OR`
#' @param upstream_steps non-negative integer. If `genes` is specified, include
#' genes up to this many directed upstream steps before subsetting
#' @import OmnipathR
#' @importFrom dplyr %>%
#' @importFrom igraph delete.edges feedback_arc_set
#' @import tidygraph
#' @examples
#' extractDAGfromOmnipath("TGF")
#' @export
#' @return list of igraph object and omitted arcs
extractDAGfromOmnipath <- function(genes=NULL, organism="mouse",
    pathway_name=NULL, target_source="AND", mfas=TRUE, upstream_steps=0) {
    if (is.null(pathway_name) & is.null(genes)) {
        stop("Please specify pathway name or gene set")
    }
    if (!is.numeric(upstream_steps) || length(upstream_steps) != 1 ||
        is.na(upstream_steps) || upstream_steps < 0 ||
        upstream_steps != as.integer(upstream_steps)) {
        stop("upstream_steps should be a single non-negative integer")
    }
    upstream_steps <- as.integer(upstream_steps)

    doro <- dorothea(organism=organism)
    
    ## Pathway subsetting
    if (!is.null(pathway_name)) {
        network_slk_pw <- annotated_network(doro, 'SignaLink_pathway')
        pathway_list <- unique(network_slk_pw$pathway_target)
        pathway_list <- pathway_list[!is.na(pathway_list)]
        if (!(pathway_name %in% pathway_list)) {
            stop("No SignaLink pathway of that name")
        }
        if (target_source=="OR") {
            sub <- network_slk_pw %>%
            filter(pathway_target==pathway_name | pathway_source==pathway_name) %>%
            filter(dorothea_level %in% c("A"))
        } else if (target_source == "AND") {
            sub <- network_slk_pw %>%
                filter(pathway_target==pathway_name & pathway_source==pathway_name) %>%
                filter(dorothea_level %in% c("A"))       
        } else {
            stop("target_source should be OR or AND")
        }
    } else {
        sub <- doro
    }
    
    ## Gene subsetting
    if (!is.null(genes)) {
        genes <- .expand_with_upstream_genes(sub, genes, upstream_steps)
        if (target_source=="OR") {
            sub <- sub %>% 
            filter(source_genesymbol %in% genes | target_genesymbol %in% genes)            
        } else if (target_source == "AND") {
            sub <- sub %>% 
            filter(source_genesymbol %in% genes & target_genesymbol %in% genes)            
        } else {
            stop("target_source should be OR or AND")
        }
    }
    
    if (dim(sub)[1]==0) {stop("No interaction remained.")}
    
    graph <- sub[,c("source_genesymbol","target_genesymbol")] %>%
        `colnames<-`(c("from","to"))
    
    graph.nodup <- graph[!duplicated(graph[, c("from","to")]),]
    ig <- igraph::graph_from_data_frame(graph.nodup)

    if (isTRUE(mfas)) {
        fas <- feedback_arc_set(ig)
        ig2 <- delete.edges(ig,  fas)
        if (!is_dag(ig2)) {stop("The resulting graph is not DAG")}
    } else {
        ig2 <- ig
    }
    
    comp <- to_components(ig2)
    maxind <- which.max(unlist(lapply(comp, function(XX) {
        dm <- XX %N>% data.frame(); dim(dm)[1]})))
    ig3 <- comp[[maxind]]
    
    if (isTRUE(mfas)) {
        return(list("igraph"=ig3, "fas"=fas))
    } else {
        return(ig3)
    }
}

.expand_with_upstream_genes <- function(sub, genes, upstream_steps) {
    genes <- unique(as.character(genes))
    if (upstream_steps == 0 || nrow(sub) == 0) {
        return(genes)
    }

    graph <- sub[, c("source_genesymbol", "target_genesymbol")]
    graph <- graph[stats::complete.cases(graph), , drop = FALSE]
    if (nrow(graph) == 0) {
        return(genes)
    }
    colnames(graph) <- c("from", "to")

    ig <- igraph::graph_from_data_frame(graph, directed = TRUE)
    matched_genes <- intersect(genes, igraph::V(ig)$name)
    if (length(matched_genes) == 0) {
        return(genes)
    }

    upstream <- igraph::neighborhood(
        ig,
        order = upstream_steps,
        nodes = matched_genes,
        mode = "in"
    )
    upstream_genes <- unique(unlist(lapply(upstream, function(node_ids) {
        igraph::V(ig)$name[node_ids]
    })))

    unique(c(genes, upstream_genes))
}

