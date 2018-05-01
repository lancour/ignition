#' Creates an adjacency matrix from a set of edges.
#' 
#' The first step in the ignition workflow (if you are using your own edge set) is 
#' to convert the edges into an adjacency matrix. This method provides a simple
#' framework for doing so.
#' 
#' @param edge.set A data frame with two columns. Each row is an interaction and each column
#' is a gene involved in the interaction. Optionally, the data frame can additionally contain
#' a third column containing a numeric confidence measure for each interaction. Interactions
#' with missing or negative confidence measures are ignored.
#' @return An adjacency matrix representation of the given edge set.
#' @examples
#' data(ignition.example.edges)
#' adj.mat = CreateAdjMatrix(ignition.example.edges)
#' @export
CreateAdjMatrix <- function(edge.set){
  weighted = F
  if (ncol(edge.set) == 3){
    weighted = T
    weight.vec = as.numeric(edge.set[,3])
    keep.ind = ((!is.na(weight.vec)) & (weight.vec > 0))
    edge.set = edge.set[keep.ind,]
    weight.vec = as.numeric(edge.set[,3])
  }
  
  e1 = as.character(toupper(edge.set[,1]))
  e2 = as.character(toupper(edge.set[,2]))
  node.set = unique(c(e1,e2))
  gene.hash = hash::hash(node.set, 1:length(node.set))
  adj.mat = matrix(0.0,length(gene.hash),length(gene.hash))
  for(i in 1:length(e1)){
    par1 = gene.hash[[e1[i]]]
    par2 = gene.hash[[e2[i]]]
    if (!is.null(par1) && !is.null(par2)){
      if (par1 != par2){
        if (weighted == T){
          adj.mat[par1,par2] = weight.vec[i]
          adj.mat[par2,par1] = weight.vec[i]
        } else{
          adj.mat[par1,par2] = 1
          adj.mat[par2,par1] = 1
        }
      }
    }
  }
  rownames(adj.mat) = node.set
  colnames(adj.mat) = node.set
  return(adj.mat)
}
