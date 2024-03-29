% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processing_wrappers.R
\name{integrate_data}
\alias{integrate_data}
\title{Integrate data from IL22+ and IL22- conditions}
\usage{
integrate_data(data = NULL)
}
\arguments{
\item{data}{The output from \code{\link{preprocess_data}} which can
be directly fed in using the pipe operator \code{\%>\%}}
}
\value{
An integrated Seurat object to be passed to \code{\link{process_data}}
}
\description{
This function takes the output from \code{\link{preprocess_data}} and
performs data integration as described in
\href{https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8}{ Stuart,
Butler et al, 2019 }.
}
\details{
Follows the steps documented in the Seurat vignette found
\href{https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1}{here}.
Briefly, \enumerate{ \item Calls \code{\link[Seurat]{SplitObject}} to
subset the object by condition and outputs to a list \item Calls the
function \code{\link[Seurat]{SCTransform}} to perform SCTransform
normalization (described in \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1}{Hafemeister and Satija, 2019}) on the split objects
\item Calls \code{\link[Seurat]{SelectIntegrationFeatures}} to select integration features (n = 3000)
\item Calls \code{\link[Seurat]{PrepSCTIntegration}} to run pre-integration checks
\item Calls \code{\link[Seurat]{FindIntegrationAnchors}} to define integration anchors using the integration features selected above
\item Calls \code{\link[Seurat]{IntegrateData}} to perform the data final data integration}
}
\section{Note}{
 Should you receive the error below, you need to increase the
size limit for global variables accessed by functions passed in the
\code{future} package. This can be done by setting
\code{options(future.globals.maxSize = 1000 * 1024^2)}, for example (1GB
limit). Be aware that this will increase RAM usage, so use carefully. More
information can be found
\href{https://satijalab.org/seurat/articles/future_vignette.html}{here}.

\code{The total size of the X globals that need to be exported for the
  future expression ('FUN()') is X GiB. This exceeds the maximum allowed size
  of 500.00 MiB (option 'future.globals.maxSize'). The X largest globals are
  ... }
}

\examples{
\dontrun{
preprocess_data() \%>\% integrate_data()}
}
\seealso{
Other data processing: 
\code{\link{preprocess_data}()},
\code{\link{process_data}()}
}
\concept{data processing}
