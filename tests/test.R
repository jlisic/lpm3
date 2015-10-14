library(lpm3)
library(BalancedSampling)

#set.seed(100)

N <- 100000
n <- 100 
x <- cbind( runif(N), runif(N)) 



Cprog <- proc.time()
set.seed(100)
sampled2 <- lpm2( rep(n/N,N),x  )
print("lpm2 running time")
print(proc.time() - Cprog) 

Cprog <- proc.time()
set.seed(100)
sampled3 <- lpm3( rep(n/N,N),x  )
print("lpm3 running time")
print(proc.time() - Cprog) 

print(max( sampled2 - sampled3))

# iteration 1
#set.seed(100)
#r1 <- runif(N)
#r2 <- runif(N)
#
#i <- 0 
#
#j <- floor(r1[1] * (N - i))
#J <- j + 1
#
#x0 <- x[J,]
#w<- colSums((t(x) - x0)^2)
#w[J] <- Inf
#k <- which.min(w)
#W <- w[k] 
#cat(sprintf("k = %d, dist = %f\n",k-1,W))
#


if( F ) {
par(mfrow=c(1,2))
plot(x, cex=.1, col='orange')
points( x[sampled2,], col='black', cex=.3) 

plot(x, cex=.1, col='orange')
points( x[sampled3,], col='black', cex=.3) 
}

# before iter
#0:      0       1       2       3       4       5       6       7       8       9
#0:      0.5000  0.5000  0.5000  0.5000  0.5000  0.5000  0.5000  0.5000  0.5000  0.5000

#1:      3       1       2       0       4       5       6       7       8       9
#1:      0.5000  0.5000  0.5000  0.0000  0.5000  0.5000  0.5000  0.5000  0.5000  1.0000

#2:      3       9       2       0       4       5       6       7       8       1
#2:      0.5000  0.5000  0.5000  0.0000  0.5000  0.5000  0.5000  0.5000  0.5000  1.0000

#3:      3       9       6       0       4       5       2       7       8       1
#3:      0.5000  0.5000  0.0000  0.0000  0.5000  0.5000  1.0000  0.5000  0.5000  1.0000

#4:      3       9       6       0       4       5       2       7       8       1
#4:      1.0000  0.5000  0.0000  0.0000  0.5000  0.0000  1.0000  0.5000  0.5000  1.0000

#5:      3       9       6       0       2       5       4       7       8       1
#5:      1.0000  0.5000  0.0000  0.0000  0.5000  0.0000  1.0000  0.5000  0.5000  1.0000

#6:      3       9       6       0       2       7       4       5       8       1
#6:      1.0000  0.5000  0.0000  0.0000  0.5000  0.0000  1.0000  0.0000  1.0000  1.0000

#7:      3       9       6       0       2       7       1       5       8       4
#7:      1.0000  1.0000  0.0000  0.0000  0.0000  0.0000  1.0000  0.0000  1.0000  1.0000

#8:      3       9       6       0       2       7       1       8       5       4

#8:      1.0000  1.0000  0.0000  0.0000  0.0000  0.0000  1.0000  0.0000  1.0000  1.0000



