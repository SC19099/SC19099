#' @title Generate the Gram matrix using R
#' @description Generate the Gram matrix which consists of the inner products of vectors
#' @param m Input matrix
#' @return a matrix consisting of inner products
#' @examples
#' \dontrun{
#' data(data)
#' g <- Gram(state.x77)
#' }
#' @export
Gram <- function(m){
  n <- nrow(m)
  G <- matrix(nrow = n,ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      G[i,j] <- m[i,] %*% m[j,] 
    }
  }
  return(G)
}