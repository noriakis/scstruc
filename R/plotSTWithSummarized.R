#' plotSTWithSummarized
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
plotSTWithSummarized <- function(spe, gsv, from_node, to_node, label="label",
	sample_id=NULL, image_id=NULL) {
	stopifnot(
        is(spe, "SpatialExperiment")
	)
	        
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

	img_raster <- imgRaster(spe, sample_id, image_id)
	dim_raster <- img_raster |> dim()

	## Obtain raw coordinates
	raw <- spatialCoords(spe) |> data.frame()
	raw <- cbind(raw, colData(spe)[,c("in_tissue")],
		colData(spe)[[label]])
	raw <- raw |> `colnames<-`(c("x","y","in_tissue","label"))
	raw$x <- raw$x * sf
	raw$y <- raw$y * sf
	raw <- data.frame(raw)
	raw$y <- -1 * raw$y
	raw$in_tissue <- as.logical(raw$in_tissue)
	
	plot_num <- gsv |> filter(.data$from==from_node) |> 
	    filter(.data$to==to_node) |>
	    select(ends_with("summarize"))
	plot_num <- plot_num[1,] |> as.numeric() |>
    	setNames(gsub("_summarize","",colnames(plot_num)))
    	
    
    raw <- raw |> mutate(value=plot_num[label])
    
    ## Return plot
    ggplot(raw, layout="manual", x=x, y=y) +
        annotation_raster(img_raster,
            xmin=0, xmax=dim_raster[1],
            ymin=-1*dim_raster[2], ymax=0)+
    geom_point(aes(x=x, y=y, fill=value),
               alpha=0.8, shape=21, size=1, color="grey20")+
    coord_fixed(xlim=c(0,dim_raster[1]), ylim=c(-1*dim_raster[2],0))+
    ggtitle(paste0(from_node,"->",to_node))+
    scale_fill_gradient(low="white", high="red")
}