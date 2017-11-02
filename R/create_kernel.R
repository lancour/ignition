#' Computes the regularized laplacian matrix from an adjacency matrix
#' 
#' This method takes an adjacency matrix, which can be created with 
#' \code{\link{CreateAdjMatrix}}, and computes the diffusion kernel matrix. If the
#' dimensions of the adjacency matrix are large (at least 10,000 x 10,000), then this
#' step can take an hour or more to run.
#'  
#' @param adj.mat Required. An adjacency matrix. Can be produced using \code{\link{CreateAdjMatrix}}
#' @param lambda Optional. Defaults to standard of 0.1. Adjusts the amount of diffusion done. 
#' Not recommended to change unless there is specific rationale for doing so.
#' @param autosave Optional. Since this function can take a while to compute, it may be
#' preferable to have the kernel be saved immediately once it is computed. If autosave is set
#' to TRUE, a copy of the kernel will be saved in the current working directory.
#' @return The regularized laplacian kernel matrix.
#' @examples
#' data(ignition.example.edges)
#' adj.matrix = CreateAdjMatrix(ignition.example.edges)
#' kernel = CreateKernel(adj.matrix)   #if not using autosave
#' kernel = CreateKernel(adj.matrix,autosave=TRUE)    #if using autosave
#' @export
CreateKernel <- function(adj.mat,lambda=0.1,autosave=FALSE){
  laplace.mat = -1 * abs(adj.mat)
  for (i in 1:nrow(laplace.mat)){
    laplace.mat[i,i] = 0
    laplace.mat[i,i] = -1 * sum(laplace.mat[i,])
  }
  id.mat = diag(nrow(laplace.mat))
  adjust.mat = lambda * laplace.mat
  rm(laplace.mat)
  adjust.mat = adjust.mat + id.mat
  rm(id.mat)
  chol.mat = chol(adjust.mat)
  rm(adjust.mat)
  ker.mat = chol2inv(chol.mat)
  rm(chol.mat)
  rownames(ker.mat) = rownames(adj.mat)
  colnames(ker.mat) = colnames(adj.mat)
  if (autosave==TRUE){
    saveRDS(ker.mat,paste0("kernel_",lambda,"_thresh.rds"))
  }
  return(ker.mat)
}
