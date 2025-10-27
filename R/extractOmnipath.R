#' @title extractDAGfromOmnipath
#' This function takes the pathway name in SignaLink and
#' 1. Subset the pathway-related genes in dorothea interaction and remove duplicates
#' 2. Identify and remove minimal feedback arc set (make DAG)
#' 3. Takes the largest components in the graph
#' @param mfas force DAG
#' @import OmnipathR
#' @importFrom dplyr %>%
#' @importFrom igraph delete.edges feedback_arc_set
#' @import tidygraph
#' @examples
#' extractDAGfromOmnipath("TGF")
#' @export
extractDAGfromOmnipath <- function(genes=NULL,
    pathway_name=NULL, target_source="AND", mfas=TRUE) {
    if (is.null(pathway_name) & is.null(genes)) {
        stop("Please specify pathway name or gene set")
    }

    doro <- dorothea()
    
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

