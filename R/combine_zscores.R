#' Combines two Z scores for each gene
#' 
#' If a user has some form of gene-level significance data, for example from a genetic case
#' control analysis, they can use this method to combine those measures with network Zscores
#' 
#' @param gene.names Required. Gene identifiers for each row of the inputted Zscores. Should
#' be the same length as gene.z and net.z.
#' @param gene.z Required. Z scores for each gene taken from significance in a case control
#' genetic study, for example. Note that these should be one directional Z scores, such that 
#' positive Z scores are highly significant and negative Z scores are highly insignificant.
#' @param net.z Required. Z scores for each gene produced by \code{\link{GeneratePredictions}}
#' @param gene.weight Optional. Defaults to 0.73 as determined in our original paper. This is
#' the weight that is assigned to gene.z scores (and then net.z scores is assigned a weight 
#' of (1 - gene.weight)
#' @return A data frame containing the combined Z scores for each gene.
#' @examples
#' data(ignition.example.edges)
#' data(ignition.example.genetic)
#' known.gene.set = c('B', 'I')
#' adj.matrix = CreateAdjMatrix(ignition.example.edges)
#' kernel = CreateKernel(adj.matrix)
#' net.predictions = GeneratePredictions(kernel, known.gene.set)
#' ignition.example.genetic$gwasz = qnorm(ignition.example.genetic$pvals, lower.tail = FALSE)
#' merged.data.frame = merge(net.predictions, ignition.example.genetic, by = "gene")
#' CombineZscores(merged.data.frame$gene,merged.data.frame$gwasz,merged.data.frame$netz)
#' @export
CombineZscores <- function(gene.names, gene.z, net.z, gene.weight=0.73){
  net.weight = 1.0 - gene.weight
  new.z = (gene.z * gene.weight) + (net.z * net.weight)
  norm.denom = ((gene.weight ^ 2) + (net.weight ^ 2)) ^ (1/2)
  new.z = new.z / norm.denom
  p.vals = pnorm(new.z,lower.tail=F)
  ret.df = data.frame("gene" = gene.names, "z" = new.z, "p" = p.vals)
  return(ret.df)
}
