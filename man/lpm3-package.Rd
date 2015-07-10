\name{saAlloc-package}
\alias{saAlloc-package}
\docType{package}
\title{
A Fast K-D Tree Implementation of the Local Pivotal Method 
}
\description{
This program performs the local pivotal method for balanced sampling using k-d trees obtaining O( log(n) n ) computatinal complexity. 
}
\details{
\tabular{ll}{
Package: \tab lpm3\cr
Type: \tab Package\cr
Version: \tab 0.5\cr
Date: \tab 2015-07-10\cr
License: \tab BSD 2-clause License\cr
}

lpm3: currently does not do correct handling of ties on searches.

}
\author{
Jonathan Lisic

Maintainer: <jonathan.lisic@nass.usda.gov>
}
\keyword{ package }
\keyword{ sampling survey }
\examples{
# balanced sample with sampling rate at 0.10

N <- 1000
n <- 100
x <- cbind( runif(N), runif(N)) 

Cprog <- proc.time()
sampled3 <- lpm3( rep(n/N,N),x  )
print("lpm3 running time")
print(proc.time() - Cprog) 

print(length(sampled3))

}
