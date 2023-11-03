#' plotSTNet
#'
#' Overlay igraph network on histology image
#' 
#' @param spe SpatialExperiment
#' @param ig igraph or data.frame object containing `from` and `to` column
#' @param label label in `colData(spe)` that matches the node labels in ig
#' @param edge_width edge width in igraph or data.frame
#' @param edge_color edge color in igraph or data.frame
#' @param edge_range edge width range
#' @param sample_id sample_id in `spe`, used to subset image
#' @param image_id image_id in `spe`, used to subset image
#' @param return_tbl_graph return tbl graph only
#' @param point_label background point label
#' @param cell_cluster_label label of each node, should be named vector
#' @param directed network is directed or not
#' @param use_ggfx use ggfx or not
#' @param point_size representative point size
#' @param rep_colour representative point colour
#' @param edge_base_colour fixed edge colour
#' @param edge_base_width fixed edge width
#' @export
#' @return ggplot
#' 
plotSTNet <- function(spe, ig, label="label", edge_width=NULL, edge_color=NULL,
	edge_range=c(0.5,1.5), sample_id=NULL, image_id=NULL, return_tbl_graph=FALSE,
	point_label="ground_truth", cell_cluster_label=NULL, directed=TRUE,
	use_ggfx=TRUE, point_size=2, rep_colour="black", edge_base_colour="grey50",
	edge_base_width=0.5) {
	
	stopifnot(
        is(spe, "SpatialExperiment")
	)
	
	if (is.data.frame(ig)) {
		if (!"from" %in% colnames(ig)) stop("No from column in data frame") 
		if (!"to" %in% colnames(ig)) stop("No to column in data frame") 
	} else {
        if (!"igraph" %in% class(ig)) {stop("please provide igraph object")}	
	}
    coords <- spatialCoords(spe)
    
    ## Check duplicate rows
    coords <- coords[!duplicated(coords), ] |> data.frame() |> `colnames<-`(c("x","y"))
    
    if (!is.null(label)) {
    	coords <- cbind(coords, colData(spe)[["label"]]) |>
    	    `colnames<-`(c("x","y","label"))    	
    } else {
    	coords$label <- row.names(coords)
    }
        
    ## scale factor based on ID
    if (!is.null(sample_id) & !is.null(image_id)) {
	    sf <- imgData(spe) |> 
	        data.frame() |>
	        dplyr::select(c("sample_id","image_id","scaleFactor")) |>
	        dplyr::filter(.data$sample_id==sample_id) |>
	        dplyr::filter(.data$image_id==image_id) |>
	        pull(scaleFactor)
	} else {
		sf <- imgData(spe)[1,]$scaleFactor
	}

	nodes <- coords |>
	    group_by(label) |>
	    summarise(x=mean(x), y=mean(y)) |> 
	    mutate(name=label) |>
	    mutate(y=-1*y) |> 
	    mutate(x=x*sf, y=y*sf)
	if (!is.null(cell_cluster_label)) {
		nodes <- nodes |>
		    mutate(cell_cluster_label=cell_cluster_label[name])
	}
	
	if (!is.data.frame(ig)) {
    	elist <- ig |> igraph::as_data_frame()
	}
	plotG <- tbl_graph(nodes=nodes, edges=ig)
	if (return_tbl_graph) {
        return(plotG)
	}

	img_raster <- imgRaster(spe, sample_id, image_id)
	dim_raster <- img_raster |> dim()

	## Obtain raw coordinates
	raw <- spatialCoords(spe) |> data.frame()
	if (is.null(point_label)) {
		raw <- cbind(raw, colData(spe)[,c("in_tissue")]) |>
		    data.frame() |> `colnames<-`(c("x","y","in_tissue"))
	} else {
		raw <- cbind(raw, colData(spe)[,c("in_tissue", point_label)]) |>
		    data.frame() |> `colnames<-`(c("x","y","in_tissue", point_label))
	}
	raw$x <- raw$x * sf
	raw$y <- raw$y * sf
	raw <- data.frame(raw)
	raw$y <- -1 * raw$y
	raw$in_tissue <- as.logical(raw$in_tissue)
    
    
    ## Ensure the variables
    if (directed) {
    	arrow <- arrow(type="closed", length=unit(1,"mm"))
    } else {
    	arrow <- NULL
    }
    if (is.null(point_label)) {
    	point_fill <- NULL
    } else {
    	point_fill <- quo(.data[[point_label]])
    }
    if (!is.null(edge_width)) {
    	.edge_width <- quo(.data[[edge_width]])
    } else {
    	.edge_width <- NULL
    }
    if (!is.null(edge_color)) {
    	.edge_color <- quo(.data[[edge_color]])
    } else {
    	.edge_color <- NULL
    }
    .edge_fnc <- function(use_ggfx=FALSE, geom="geom_edge_link") {
    	args <- list()
    	args[["mapping"]] <- aes(width=!!.edge_width, color=!!.edge_color)
    	args[["arrow"]] <- arrow
    	args[["end_cap"]] <- circle(3,"mm")
    	args[["start_cap"]] <- circle(3,"mm")
    	if (is.null(edge_color)) {
    		args[["edge_colour"]] <- edge_base_colour
    	}
    	if (is.null(edge_width)) {
    		args[["width"]] <- edge_base_width
    	}

    	if (use_ggfx) {
		    ggfx::with_outer_glow(
		        do.call(geom, args),
		        colour="white"
		    )    		
    	} else {
	        do.call(geom, args)
    	}
    }
    
    ## Return plot
    ggraph(plotG, layout="manual", x=x, y=y) +
        annotation_raster(img_raster,
            xmin=0, xmax=dim_raster[1],
            ymin=-1*dim_raster[2], ymax=0)+
    geom_point(aes(x=x, y=y, fill=!!point_fill),
               alpha=0.5, shape=21, color="grey",
               data=raw[raw$in_tissue,])+
    coord_fixed(xlim=c(0,dim_raster[1]), ylim=c(-1*dim_raster[2],0))+
    .edge_fnc(use_ggfx)+
    scale_edge_width(range=edge_range)+
    geom_node_point(aes(fill=cell_cluster_label), size=point_size,
    	colour=rep_colour, shape=21)
}