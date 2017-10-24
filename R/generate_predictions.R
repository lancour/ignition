#' Generates network scores for all genes in a network
#' 
#' This method takes a gene to gene distance matrix and computes the summed proximity
#' of each gene to a set of genes. The set of genes is presumably a collection of known
#' genes for a phenotype of interest.
#' 
#' @param kernel Required. Produced by \code{\link{CreateKernel}} or user may provide own
#' gene to gene distance kernel. Row and column names must be set to gene symbols.
#' @param pheno.gene.set Required. A vector of previously known phenotype-related genes. WIll
#' automatically be coerced to uppercase unless overriden (see autocaps param).
#' @param autocaps Optional. Defaults to TRUE. Converts all entries in pheno.gene.set to
#' uppercase to ensure accurate mapping to genes in kernel.
#' @param pseudo.add Optional. Defaults to 0.1. Adds a small pseudocount to prevent any gene
#' from having a percentile of exactly 0 or 1 (qnorm produces an infinite Z score otherwise)
#' @return A data frame containing the computed percentiles and Z scores for each gene.
#' @examples
#' data(ignition.example.edges)
#' adj.matrix = CreateAdjMatrix(ignition.example.edges)
#' kernel = CreateKernel(adj.matrix)
#' known.gene.set = c('B', 'I')
#' GeneratePredictions(kernel, known.gene.set)
#' @export
GeneratePredictions <- function(kernel, pheno.gene.set, autocaps = FALSE, pseudo.add = 0.1){
  k.genes = rownames(kernel)
  base.vec = (k.genes %in% pheno.gene.set) * 1
  base.ind = which(base.vec ==1)
  diff.scores = kernel %*% base.vec
  print.df = data.frame("gene" = k.genes, "score" = diff.scores, stringsAsFactors=FALSE)
  print.df = print.df[-c(base.ind),]
  print.df$ptile = rank(print.df$score) / (nrow(print.df) + pseudo.add)
  print.df$netz = qnorm(1 - print.df$ptile,lower.tail=FALSE)
  print.df = print.df[order(-print.df$netz,print.df$gene),]
  return(print.df)
}
