#' Raw input data
#'
#' A list containing the raw count output data from CellRanger for each of the 2 culture conditions (IL22+ and IL22-)
#' which were read into R using the \code{Seurat::Read10X(..., strip.suffix=T)} function call
#'
#' @format A list containing 2 sparse matrices (of type 'dgCMatrix') where rows are genes and cols are cells
#' \describe{
#'   \item{IL22+}{dgCMatrix with 33538 rows (features) and 1427 columns (cells)}
#'   \item{IL22-}{dgCMatrix with 33538 rows (features) and 2168 columns (cells)}
#' }
#' @concept data
"raw_data"

#' Stressed cells
#'
#' A character vector containing the IDs of cells to be excluded based
#' on their expression of stress-response genes
#'
#' @format A character vector with 37 entries
#'
#' @concept data
#'
"stressed_cells"

#' Cell cycle genes
#'
#' A character vector containing cell-cycle phase-specific gene symbols; 1:43
#' are S-Phase specific, 44:97 are G2/M specific. Used for cell cycle scoring by
#' \code{\link{process_data}}
#'
#' @format A character vector with 97 entries
#' @concept data
#' @source \url{https://satijalab.org/seurat/articles/cell_cycle_vignette.html}
"cc_genes"

#' Raw velocity data
#'
#' Seurat object generated from the output loom files produced by the velocyto
#' command line tool. A detailed overview of how this file was generated can be
#' found in the package vignette.
#'
#' @format A Seurat object with 4 assay slots: spliced, unspliced, ambiguous and RNA
#'
#' @concept data
"velocity_raw"
