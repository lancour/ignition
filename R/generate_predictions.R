#' Generates network scores for all genes in a network
#' 
#' This method takes a gene to gene distance matrix and computes the summed proximity
#' of each gene to a set of genes. The set of genes is presumably a collection of known
#' genes for a phenotype of interest.
#' 
#' @param kernel Required. Produced by \code{\link{CreateKernel}} or user may provide own
#' gene to gene distance kernel. Row and column names must be set to gene symbols.
#' @param pheno.gene.set Required. A vector of previously known phenotype-related genes. WIll
#' automatically be coerced to uppercase unless overriden (see autocaps param).
#' @param autocaps Optional. Defaults to TRUE. Converts all entries in pheno.gene.set to
#' uppercase to ensure accurate mapping to genes in kernel.
#' @param invert.distance Optional. Defaults to FALSE. If using a different kernel other than
#' the difussion kernel produced by \code{\link{CreateKernel}}, this may need to be toggled
#' to TRUE. For example, high values in a diffusion kernel are good, whereas for other distance
#' measures, such as shortest path, low values are good.
#' @param ignore.neighbors Optional. Defaults to FALSE. If set to TRUE, then predictions will
#' only be generated for genes that do not directly interact with the pheno.gene.set. This allows
#' for investigation of indirect interactions that can be detected with ignition's diffusion
#' method. Note that if this is set to TRUE, an adjacency matrix such as that produced by
#' \code{\link{CreateAdjMatrix}} must be provided
#' @param adj.matrix Optional. An adjacency matrix of the same dimensions as the kernel. If
#' using a custom kernel, make sure that the ith row/col of the kernel correspond to the ith
#' row/col of the adjacency matrix. If the user has specified ignore.neighbors = TRUE, then 
#' this argument is required.
#' @return A data frame containing the computed percentiles and Z scores for each gene.
#' @examples
#' data(ignition.example.edges)
#' adj.matrix = CreateAdjMatrix(ignition.example.edges)
#' kernel = CreateKernel(adj.matrix)
#' known.gene.set = c('B', 'I')
#' GeneratePredictions(kernel, known.gene.set)
#' @export
GeneratePredictions <- function(kernel, pheno.gene.set, autocaps = FALSE,
                                invert.distance = FALSE, ignore.neighbors = FALSE, 
                                adj.matrix = NULL){
  if (ignore.neighbors == TRUE){
    if (!(class(adj.matrix) == "matrix")) {
      stop("If using ignore.neighbors option, a matrix object must be provided to the
           adj.matrix argument \n")
    }
  }
  k.genes = rownames(kernel)
  base.vec = (k.genes %in% pheno.gene.set) * 1
  base.ind = which(base.vec ==1)
  if (length(base.ind) > 0){
    print(paste(length(base.ind), "of", length(pheno.gene.set), "genes present in kernel"))
  } else{
    stop("None of the genes in pheno.gene.set were able to be mapped to the kernel. Please 
    make sure you are using correct Human Gene Nomenclature Symbols, and that you are
    not using old versions of the symbols. Also make sure that the rows and columns of the
    kernel are labeled correctly with HGNC symbols.")
  }
  diff.scores = kernel %*% base.vec
  print.df = data.frame("gene" = k.genes, "score" = diff.scores, stringsAsFactors=FALSE)
  
  drop.ind = base.ind
  if (ignore.neighbors == TRUE){
    for (g in 1:length(base.ind)){
      cur.ind = base.ind[g]
      drop.ind = c(drop.ind, which(adj.matrix[cur.ind, ] == 1))
    }
    drop.ind = unique(drop.ind)
  }
  
  print.df = print.df[-c(drop.ind),]
  print.df$ptile = rank(print.df$score) / nrow(print.df)
  adjust.val = min(print.df$ptile) / 2.0
  print.df$ptile = print.df$ptile - adjust.val
  if (invert.distance == TRUE){
    print.df$netz = qnorm(print.df$ptile,lower.tail=FALSE)
  } else{
    print.df$netz = qnorm(1 - print.df$ptile,lower.tail=FALSE)
  }
  print.df = print.df[order(-print.df$netz,print.df$gene),]
  return(print.df)
}
