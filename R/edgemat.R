#' @title Generate the edge matrix using R
#' @description Generate relationships between variables with a set probability
#' @param p Variable dimension
#' @param pr Probability matrix
#' @return a matrix consisting of edge information
#' @examples
#' \dontrun{
#' p <- 100
#' x <- runif(p,0,1)
#' y <- runif(p,0,1)
#' node <- cbind(x,y)
#' distance <- as.matrix(dist(node))
#' pr <- matrix(0,nrow = p,ncol = p)
#' for (i in 1:(p-1)) {
#' for (j in (i+1):p) {
#' d <- distance[i,j]
#' pr[i,j] <- dnorm(d * sqrt(p))
#' }}
#' e <- edgemat(p,pr)
#' @export
edgemat <- function(p,pr){
  mat <- matrix(0,nrow = p,ncol = p) 
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      prob <- pr[i,j]
      lineornot <- sample(c(0,1),1,prob = c(1-prob,prob))
      if(lineornot == 1){
        mat[i,j] <- 1
      }
    }
  }
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
}