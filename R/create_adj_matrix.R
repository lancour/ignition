#' @export
CreateAdjMatrix <- function(edge.set){
  e1 = as.character(toupper(edge.set[,1]))
  e2 = as.character(toupper(edge.set[,2]))
  node.set = unique(c(e1,e2))
  gene.hash = hash::hash(node.set, 1:length(node.set))
  adj.mat = matrix(0,length(gene.hash),length(gene.hash))
  for(i in 1:length(e1)){
    par1 = gene.hash[[e1[i]]]
    par2 = gene.hash[[e2[i]]]
    if (!is.null(par1) && !is.null(par2)){
      if (par1 != par2){
        adj.mat[par1,par2] = 1
        adj.mat[par2,par1] = 1
      }
    }
  }
  rownames(adj.mat) = node.set
  colnames(adj.mat) = node.set
  return(adj.mat)
}
