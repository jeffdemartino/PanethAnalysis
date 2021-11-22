#' Raw input data
#'
#' A list containing the raw count data from CellRanger for each of the 2 culture conditons (IL22+ and IL22-)
#' which were read into R using the \code{Seurat::Read10X(..., strip.suffix=T)} function call
#'
#' @format A list containing 2 sparse matrices (of type 'dgCMatrix') where rows are genes and cols are cells
#' \describe{
#'   \item{IL22+}
#'   \item{IL22-}
#'   ...
#' }
#' @concept data
#'
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
#'
"stressed_cells"
