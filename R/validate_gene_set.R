#' @export
ValidateGeneSet <- function(kernel, pheno.gene.set){
  ker.genes = rownames(kernel)
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