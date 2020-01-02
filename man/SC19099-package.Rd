\name{SC19099-package}
\alias{SC19099-package}
\alias{SC19099}
\docType{package}
\title{
  Compute the inner product matrix and edge matrix
}
\description{
  Compute the inner product matrix of a set of vectors and edge matrix which shows the relationship of variables.
}
\details{
  SC19099 is a simple R package developed to implement two functions. Two functions are _Gram_ (generate a matrix consisting of the inner product of vectors) and _edgemat_ (establish a relationship matrix between variables with a set probability) respectively. Specifically, by inputting a set of vectors, the internal product between the vectors can be calculated using _Gram_, and these internal products form a symmetric matrix; By entering the set probability and the dimension of the variable, the probability can be used to determine whether the variables are related by _edgemat_, If relevant, matrix corresponding element is 1, otherwise 0.
}
\author{
Yu Xia
Maintainer: Yu Xia <xiayu01@mail.ustc.edu.cn>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
     ## Optional simple examples of the most important functions
     ## These can be in \dontrun{} and \donttest{} blocks.   
  }
}
