% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing_wrappers.R
\name{process_data}
\alias{process_data}
\title{Process integrated data}
\usage{
process_data(data = NULL, output.dir = ".")
}
\arguments{
\item{data}{The output from \code{\link{integrate_data}} which can
be directly fed in using the pipe operator \code{\%>\%}}

\item{output.dir}{Path to directory where DGE output should be saved
(default is ".")}
}
\value{
An integrated Seurat object for which PCA and UMAP have been computed as well as (annotated) graph-based clustering
}
\description{
Runs dimensional reduction, graph-based clustering and cluster annotation
}
\details{
This function runs the following data processing steps:
\enumerate{\item Performs PCA by calling \code{\link[Seurat]{RunPCA}}
\item Runs UMAP using the first 25 principal components (see \code{\link[Seurat]{RunUMAP}})
\item Performs graph-based clustering by calling \code{\link[Seurat]{FindNeighbors}} (using the first 25 principal components) and \code{\link[Seurat]{FindClusters}} (resolution of 0.8, Leiden algorithm)
\item Calls \code{\link[Seurat]{FindAllMarkers}} to perform differential gene expression analysis between the clusters defined in the previous step (only.pos = T, logfc.threshold = 0.5). Results are written in the directory specified by \code{output.dir} as "DE_perCluster.txt" (Supplemental Table 2)
\item Runs \code{\link[Seurat]{CellCycleScoring}} to infer cell cycle phase of each cell (relies on genes stored in \code{\link{cc_genes}})
\item Score cells for the Paneth cell marker gene module (DEFA5, DEFA6, PLA2G2A, PRSS2, REG3A and ITLN2)
\item Renames clusters according to manual annotations
\item Sub-clusters the Secretory progenitor cluster by calling \code{\link[Seurat]{FindSubCluster}} (resolution 1, Leiden algorithm)
\item Call Paneth cells (Paneth score > 1.5) and Tuft cells (sub-cluster 6)}
}
\examples{
\dontrun{
preprocess_data() \%>\% integrate_data() \%>\% process_data()}
}
\seealso{
Other data processing: 
\code{\link{integrate_data}()},
\code{\link{preprocess_data}()}
}
\concept{data processing}
