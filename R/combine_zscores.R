#' @export
CombineZscores <- function(gene.z,net.z,gene.weight=0.27){
  net.weight = 1.0 - gene.weight
  new.z = (gene.z * gene.weight) + (net.z * net.weight)
  norm.denom = ((gene.weight ^ 2) + (net.weight ^ 2)) ^ (1/2)
  new.z = new.z / norm.denom
  p.vals = pnorm(new.z,lower.tail=F)
  ret.df = data.frame("z" = new.z, "p" = p.vals)
  return(ret.df)
}
