#' @export
create_adj_matrix <- function(edge_set){
    e1 = as.character(toupper(edge_set[,1]))
    e2 = as.character(toupper(edge_set[,2]))
    node_set = unique(c(e1,e2))
    gene_hash = hash::hash(node_set, 1:length(node_set))
    adj_mat = matrix(0,length(gene_hash),length(gene_hash))
    for(i in 1:length(e1)){
        par1 = gene_hash[[e1[i]]]
        par2 = gene_hash[[e2[i]]]
        if (!is.null(par1) && !is.null(par2)){
            if (par1 != par2){
                adj_mat[par1,par2] = 1
                adj_mat[par2,par1] = 1
            }
        }
    }
    rownames(adj_mat) = node_set
    colnames(adj_mat) = node_set
    return(adj_mat)
}
