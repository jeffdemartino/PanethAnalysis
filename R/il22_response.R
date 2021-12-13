#' Run the IL22-response analyses
#'
#' Reproduces the IL22 response analyses used to generate XXX
#'
#' @param data The final Seurat object produced by running the
#'   \code{\link{process_data}} function
#' @param output_dir Path to directory where output files should be saved
#'   (default is "./")
#' @details This function first runs differential expression analysis between
#'   IL22+ and IL22- cells by calling the \code{\link[Seurat]{FindMarkers}}
#'   function (Wilcoxon Rank Sum test), and saves the output. The top 50 genes
#'   (by Avg Log2FC) up-regulated in IL22+ cells are then used to calculate an
#'   IL22-response module score. Finally, the differential gene expression
#'   analysis between IL22+ and IL22- cells is repeated independently for each
#'   cluster and the output is saved.
#'
#' @return This function returns the input Seurat object with IL22 response
#'   score added to the metadata \code{col.name = "IL22_Module"}. Additionally,
#'   2 output files are saved to the specified output directory: \itemize{\item
#'   "DE_IL22.txt" is a table containing the IL22+ vs IL22- differential
#'   expression results \item "DE_IL22_perCluster.txt" is a table containing the
#'   IL22+ vs IL22- differential expression results \strong{PER CLUSTER}}
#' @export
#' @concept analysis
#'
#' @examples
#' \dontrun{
#' process_data() %>% il22_response()}

il22_response <- function(data=NULL, output_dir="./"){
  Seurat::Idents(data) <- "condition"
  markers <- Seurat::FindMarkers(data, ident.1 = "IL22+", only.pos = F, logfc.threshold = 0)
  markers <- tibble::rownames_to_column(markers, var="gene")
  markers <- cbind.data.frame(markers, "avg_expression"=Matrix::rowMeans(data@assays$SCT@data)[markers$gene])
  write.table(markers, paste0(output_dir,"DE_IL22.txt"), quote = F)

  # Per cluster differential expression between IL22+ and IL22-
  Seurat::Idents(srat) <- "Final.IDS"

  ids <- as.character(unique(Seurat::Idents(srat)))
  markers.clus <- list()
  for (cluster in 1:length(ids)){
    srat.sub <- subset(srat, cells = Seurat::WhichCells(srat,idents = ids[cluster]))
    Seurat::Idents(srat.sub) <- "condition"
    markers.clus[[ids[cluster]]] <- Seurat::FindMarkers(srat.sub, ident.1 = "IL22+", ident.2 = "IL22-", only.pos=F,
                                                        min.cells.group = 0)
  }

  markers.clus <- lapply(markers.clus, tibble::rownames_to_column, var = "gene")

  write.table(dplyr::bind_rows(markers.clus, .id="cluster"), paste0(output_dir,"DE_IL22_perCluster.txt"), quote = F)

  # Make an IL22 induced score by taking the top 50 upregulated genes (by log2FC)
  il22.genes <- markers %>% dplyr::top_n(n = 50, wt = avg_log2FC) %>% dplyr::pull(gene)
  il22.score <- testModuleScore(Seurat::GetAssayData(data), genes.list = list(il22.genes), uniqueOnly = F, ctrl.size = 100)
  srat <- Seurat::AddMetaData(data, metadata = il22.score, col.name = "IL22_Module")

  ### IMPLEMENT IF THIS GETS ADDED TO THE MANUSCRIPT
  # Make an IL22 induced NEGATIVE score by taking the top 50 upregulated genes (by log2FC)
  #il22.genes <- markers %>% dplyr::top_n(n = -50, wt = avg_log2FC) %>% dplyr::pull(gene)
  #il22.score <- testModuleScore(Seurat::GetAssayData(data), genes.list = list(il22.genes), uniqueOnly = F, ctrl.size = 100)
  #srat <- Seurat::AddMetaData(data, metadata = il22.score, col.name = "IL22_Module_DOWN")


  return(srat)
}


