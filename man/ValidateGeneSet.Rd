% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate_gene_set.R
\name{ValidateGeneSet}
\alias{ValidateGeneSet}
\title{Conducts cross validation to test for diffusion proximity of a gene set}
\usage{
ValidateGeneSet(kernel, pheno.gene.set, autocaps = TRUE)
}
\arguments{
\item{kernel}{Required. A square, symmetrical distance kernel. \code{\link{CreateKernel}}
will produce this but the user may provide their own kernel. If you are using your own
kernel, make sure to have set the row and column names of the matrix to gene symbols.}

\item{pheno.gene.set}{Required. A character vector of genes that represent the previously
known phenotype-related genes. The function will test for proximity between these genes. 
This set will automatically be coerced to uppercase unless overriden (see autocaps param)}

\item{autocaps}{Optional. Defaults to TRUE. Converts the provided gene set to uppercase
to coincide with the upper case format used for gene identifiers.}
}
\value{
A data frame containing the cross validation percentiles for each gene.
}
\description{
Our methodology assumes that the phenotype-related genes for the phenotype of interest
actually have a significantly closer than random proximity to one another. This function
is used to test this assumption for any set of genes.
}
\examples{
data(ignition.example.edges)
adj.matrix = CreateAdjMatrix(ignition.example.edges)
kernel = CreateKernel(adj.matrix)
nonproximal.gene.set = c('B', 'I', 'F')
proximal.gene.set = c('D', 'B', 'A')
ValidateGeneSet(kernel, nonproximal.gene.set)
ValidateGeneSet(kernel, proximal.gene.set)
}

