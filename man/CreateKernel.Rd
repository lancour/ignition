% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_kernel.R
\name{CreateKernel}
\alias{CreateKernel}
\title{Computes the regularized laplacian matrix from an adjacency matrix}
\usage{
CreateKernel(adj.mat, lambda = 0.1, normalize = FALSE, autosave = FALSE)
}
\arguments{
\item{adj.mat}{Required. An adjacency matrix. Can be produced using \code{\link{CreateAdjMatrix}}}

\item{lambda}{Optional. Defaults to standard of 0.1. Adjusts the amount of diffusion done. 
Not recommended to change unless there is specific rationale for doing so.}

\item{normalize}{Optional. Defaults to false. If set to true, the normalized version 
of the laplacian is computed prior to kernel computation.}

\item{autosave}{Optional. Since this function can take a while to compute, it may be
preferable to have the kernel be saved immediately once it is computed. If autosave is set
to TRUE, a copy of the kernel will be saved in the current working directory.}
}
\value{
The regularized laplacian kernel matrix.
}
\description{
This method takes an adjacency matrix, which can be created with 
\code{\link{CreateAdjMatrix}}, and computes the diffusion kernel matrix. If the
dimensions of the adjacency matrix are large (at least 10,000 x 10,000), then this
step can take an hour or more to run.
}
\examples{
data(ignition.example.edges)
adj.matrix = CreateAdjMatrix(ignition.example.edges)
kernel = CreateKernel(adj.matrix)   #if not using autosave
kernel = CreateKernel(adj.matrix,autosave=TRUE)    #if using autosave
}

