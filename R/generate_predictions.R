#' @export
GeneratePredictions <- function(kernel, pheno.gene.set, pseudo.add = 1){
  k.genes = rownames(kernel)
  base.vec = (k.genes %in% pheno.gene.set) * 1
  base.ind = which(base.vec ==1)
  diff.scores = kernel %*% base.vec
  print.df = data.frame("gene" = k.genes, "score" = diff.scores, stringsAsFactors=FALSE)
  print.df = print.df[-c(base.ind),]
  print.df$ptile = rank(print.df$score) / (nrow(print.df) + pseudo.add)
  print.df$z = qnorm(1 - print.df$ptile,lower.tail=FALSE)
  print.df = print.df[order(-print.df$z,print.df$gene),]
  return(print.df)
}
