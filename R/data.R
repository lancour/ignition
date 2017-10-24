#' Toy gene interaction dataset for illustration of network diffusion behavior
#'
#' A dataset containing sets of fake gene interactions. Here the genes are simply letters
#' for simplicity
#'
#' @format A data frame with 12 rows and 2 variables:
#' \describe{
#'   \item{gene.one}{first gene involved in the interaction}
#'   \item{gene.two}{second gene involved in the interaction}
#'   ...
#' }
"ignition.example.edges"
#' Toy gene-level association dataset
#'
#' A dataset containing genes and p-values
#'
#' @format A data frame with 12 rows and 2 variables:
#' \describe{
#'   \item{gene}{gene symbol, in this toy example they are just letters}
#'   \item{pvals}{example p-values for each gene}
#'   ...
#' }
"ignition.example.genetic"