#' Conducts cross validation to test for diffusion proximity of a gene set
#' 
#' Our methodology assumes that the phenotype-related genes for the phenotype of interest
#' actually have a significantly closer than random proximity to one another. This function
#' is used to test this assumption for any set of genes.
#' 
#' @param kernel Required. A square, symmetrical distance kernel. \code{\link{CreateKernel}}
#' will produce this but the user may provide their own kernel. If you are using your own
#' kernel, make sure to have set the row and column names of the matrix to gene symbols.
#' @param pheno.gene.set Required. A character vector of genes that represent the previously
#' known phenotype-related genes. The function will test for proximity between these genes. 
#' This set will automatically be coerced to uppercase unless overriden (see autocaps param)
#' @param autocaps Optional. Defaults to TRUE. Converts the provided gene set to uppercase
#' to coincide with the upper case format used for gene identifiers.
#' @return A data frame containing the cross validation percentiles for each gene.
#' @examples 
#' data(ignition.example.edges)
#' adj.matrix = CreateAdjMatrix(ignition.example.edges)
#' kernel = CreateKernel(adj.matrix)
#' nonproximal.gene.set = c('B', 'I', 'F')
#' proximal.gene.set = c('D', 'B', 'A')
#' ValidateGeneSet(kernel, nonproximal.gene.set)
#' ValidateGeneSet(kernel, proximal.gene.set)
#' @export
ValidateGeneSet <- function(kernel, pheno.gene.set, autocaps = TRUE){
  ker.genes = rownames(kernel)
  if (autocaps == TRUE){
    pheno.gene.set = as.character(toupper(pheno.gene.set))
  }
  base.vec = (ker.genes %in% pheno.gene.set) * 1
  base.ind = which(base.vec==1)
  ret.df = data.frame("gene" = ker.genes[base.ind], stringsAsFactors=FALSE)
  ret.df$ptile = vector(mode="double", length=length(base.ind))
  for (g in 1:length(base.ind)){
    cur.ind = base.ind[g]
    cur.gene = ker.genes[cur.ind]
    diff.vec = base.vec
    diff.vec[cur.ind] = 0
    pheno.ind = which(diff.vec==1)
    diff.scores = kernel %*% diff.vec
    diff.genes = ker.genes
    diff.scores = diff.scores[-c(pheno.ind)]
    diff.genes = diff.genes[-c(pheno.ind)]
    new.ind = which(diff.genes == cur.gene)
    ret.df$ptile[g] = rank(diff.scores)[new.ind] / length(diff.scores)
  }
  return(ret.df)
}