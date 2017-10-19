#' @export
create_kernel <- function(adj_mat,lambda=0.1,autosave=TRUE){
    laplace_mat = -1 * (adj_mat)
    for (i in 1:nrow(laplace_mat){
        laplace_mat[i,i] = -1 * sum(laplace_mat[i,])
    id_mat = diag(nrow(laplace_mat))
    adjust_mat = lambda * laplace_mat
    adjust_mat = adjust_mat + id_mat
    ker_mat = solve(adjust_mat)
    rownames(ker_mat) = rownames(adj_mat)
    colnames(ket_mat) = colnames(adj_mat)
    if (autosave==TRUE){
        saveRDS(ker_mat,paste0("kernel_",lambda,"_thresh.rds"))
    }
}
