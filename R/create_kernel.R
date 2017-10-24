#' @export
CreateKernel <- function(adj.mat,lambda=0.1,autosave=FALSE){
  laplace.mat = -1 * (adj.mat)
  for (i in 1:nrow(laplace.mat)){
    laplace.mat[i,i] = -1 * sum(laplace.mat[i,])
  }
  id.mat = diag(nrow(laplace.mat))
  adjust.mat = lambda * laplace.mat
  adjust.mat = adjust.mat + id.mat
  ker.mat = solve(adjust.mat)
  rownames(ker.mat) = rownames(adj.mat)
  colnames(ker.mat) = colnames(adj.mat)
  if (autosave==TRUE){
    saveRDS(ker.mat,paste0("kernel_",lambda,"_thresh.rds"))
  }
  return(ker.mat)
}
