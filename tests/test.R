library(lpm3)

N <- 1000
n <- 100
x <- cbind( runif(N), runif(N)) 


#Cprog <- proc.time()
#sampled2 <- lpm2( rep(n/N,N),x  )
#print("lpm2 running time")
#print(proc.time() - Cprog) 

Cprog <- proc.time()
sampled3 <- lpm3( rep(n/N,N),x  )
print("lpm3 running time")
print(proc.time() - Cprog) 

print(length(sampled3))

