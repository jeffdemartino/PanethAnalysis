#' Prep data for RNA velocity analysis
#'
#' Prepares spliced/unspliced count data for analysis with scvelo
#'
#' @details Takes the raw spliced/unspliced count data (pre-stored as a Seurat
#'   object, see \code{\link{velocity_raw}}) and converts it to an h5ad file
#'   suitable for loading into python and analyzing with scvelo. Subsets to
#'   include only "IL22+" cells, and transfers dimensional reductions (PCA and
#'   UMAP) from the integrated, processed, Seurat object to the h5ad file.
#'   Workflow based loosely on
#'   \href{http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html}{this
#'   Seurat vignette}
#'
#'
#' @param object The final Seurat object produced by running the
#'   \code{\link{process_data}} function
#' @param output.dir Path to directory where output files should be saved
#'   (default is ".")
#'
#' @return Doesn't return an object (only saves data to disk)
#' @export
#'
#' @examples
#' \dontrun{
#' final.output <- preprocess_data() %>%
#' integrate_data() %>%
#' process_data()
#' prep_velocity(object=final.output)}
prep_velocity <- function(object=NULL, output.dir="."){

  message("Subsetting data...")
  # Subset to include only IL22+
  object.sub <- subset(object, cells = Seurat::WhichCells(object, expression = condition == "IL22+"))
  velo <- subset(PanethAnalysis::velocity_raw, cells = colnames(object.sub))


  # Transfer dimensional reductions
  velo@reductions$umap <- object.sub@reductions$umap; velo@reductions$umap@assay.used <- "RNA"
  velo@reductions$pca <- object.sub@reductions$pca; velo@reductions$pca@assay.used <- "RNA"

  velo$Final.IDS <- as.character(object.sub$Final.IDS)

  message("Creating output file...")
  # Export as h5ad
  SeuratDisk::SaveH5Seurat(velo, filename = paste0(output.dir,"//srat_velo.h5Seurat"))
  SeuratDisk::Convert(paste0(output.dir,"//srat_velo.h5Seurat"), dest = "h5ad")
  file.remove(paste0(output.dir,"//srat_velo.h5Seurat"))
}
