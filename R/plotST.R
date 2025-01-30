#' plotST
#'
#' Overlay histology images
#' 
#' @param spe SpatialExperiment
#' @param label label in `colData(spe)`
#' @param sample_id sample_id in `spe`, used to subset image
#' @param image_id image_id in `spe`, used to subset image
#' @param point_size point size
#' @export
#' @return ggplot
#' @examples
#' # library(SpatialExperiment)
#' # example(read10xVisium, echo = FALSE)
#' # plotST(spe, "in_tissue")
plotST <- function(spe, label="label",
    sample_id=NULL, image_id=NULL, point_size=1.5) {


    stopifnot(
        is(spe, "SpatialExperiment")
    )
            
    ## scale factor based on ID
    if (!is.null(sample_id) & !is.null(image_id)) {
        sf <- SpatialExperiment::imgData(spe) |> 
            data.frame() |>
            dplyr::select(c("sample_id","image_id","scaleFactor")) |>
            dplyr::filter(.data$sample_id==sample_id) |>
            dplyr::filter(.data$image_id==image_id) |>
            pull(.data$scaleFactor)
    } else {
        sf <- SpatialExperiment::imgData(spe)[1,]$scaleFactor
    }

    img_raster <- SpatialExperiment::imgRaster(spe, sample_id, image_id)
    dim_raster <- img_raster |> dim()

    ## Obtain raw coordinates
    raw <- SpatialExperiment::spatialCoords(spe) |> data.frame()
    raw <- cbind(raw, colData(spe)[,c("in_tissue")],
        colData(spe)[[label]])
    raw <- raw |> `colnames<-`(c("x","y","in_tissue","label"))
    raw$x <- raw$x * sf
    raw$y <- raw$y * sf
    raw <- data.frame(raw)
    raw$y <- -1 * raw$y
    raw$in_tissue <- as.logical(raw$in_tissue)
    
    
    ## Return plot
    ggplot(raw, layout="manual", x=.data$x, y=.data$y) +
        annotation_raster(img_raster,
            xmin=0, xmax=dim_raster[1],
            ymin=-1*dim_raster[2], ymax=0)+
    geom_point(aes(x=.data$x, y=.data$y, fill=.data$label),
               alpha=0.8, shape=21, size=point_size, color="grey20")+
    coord_fixed(xlim=c(0,dim_raster[1]), ylim=c(-1*dim_raster[2],0))
}