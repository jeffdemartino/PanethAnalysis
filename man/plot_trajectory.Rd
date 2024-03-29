% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_trajectory}
\alias{plot_trajectory}
\title{Plot trajectory figures}
\usage{
plot_trajectory(
  seurat.obj = NULL,
  slingshot.obj = NULL,
  figure = NA,
  save = FALSE,
  save.dir = "."
)
}
\arguments{
\item{seurat.obj}{The list of Seurat objects (per trajectory) generated by
\code{\link{model_diffusion}}}

\item{slingshot.obj}{The list of slingshot objects (per trajectory) generated
by \code{\link{model_slingshot}}}

\item{figure}{ID of the figure to be reproduced. Accepted values are "ED_4B"
and "ED_4C" (See Details for description of figure panels)}

\item{save}{Save the specified plot as a pdf? (default = FALSE)}

\item{save.dir}{Path to the directory where the output plot should be saved
(default = ".")}
}
\value{
The plot specified in the "figure" argument
}
\description{
Reproduces figure panels related to differentiation trajectory modelling
}
\details{
This function reproduces the following figure panels: \itemize{\item
\strong{Extended data Figure 4B}: DiffusionMaps projection of the stem cell trajectory,
colored by cell-type cluster and cell-cycle phase \item \strong{Extended data Figure
4C}: DiffusionMaps projections of the enterocyte, secretory and Paneth
cell trajectories, colored by cell-type cluster}
}
\examples{
\dontrun{
# Make each plot and save as pdf
plot_trajectory(seurat.obj = diffusion.output, slingshot.obj = slingshot.output,
                figure = "ED_4B", save = T, save.dir = ".")
plot_trajectory(seurat.obj = diffusion.output, slingshot.obj = slingshot.output,
                figure = "ED_4C", save = T, save.dir = ".")}
}
