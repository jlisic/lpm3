\name{lpm3}
\alias{lpm3}
\title{
Local Pivotal Method 
}
\description{
  The local pivotal method provides a way to perform balanced sampling.  This
  implementation uses k-d trees to replace linear searches with k-d tree queries
  moving the computational complexity of searches from O(n) to O( log(n) ).
}
\usage{

  lpm3(
     prob,
     x, 
     m=40
     ) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prob}{
    An array of length n such that the sum of prob is equal to the sample size.
  }
  \item{x}{
    A matrix of n rows and p columns, each row is assumed to be a sampling unit. 
  }
  \item{m}{
    Max leaf size used as a terminal condition for building the k-d tree.
  }
}
\author{
  Jonathan Lisic
}
\examples{


N <- 1000
n <- 100
x <- cbind( runif(N), runif(N)) 


Cprog <- proc.time()
sampled3 <- lpm3( rep(n/N,N),x  )
print("lpm3 running time")
print(proc.time() - Cprog) 

print(length(sampled3))



}

