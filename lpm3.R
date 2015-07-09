if ( !( 'loaded' %in% ls() ) )  {
  print('lpm3.so loaded')
  dyn.load('~/src/lpm3/kdtree_lpm.so')
  loaded <- T
}

  


############################################################
# x data set
# prob probability
############################################################
lpm3 <- function(
  prob,
  x,                      # data to mode seek on
  m=40
  ) {


  n <- NROW(x)
  K <- length(x) / n
  m <- min( m, n)


  # send our data to the C program
  r.result <- .C("R_lpm3",
    as.double( t(x) ),                 # 1 data we query
    as.double( prob ),                 # 2 probability
    as.integer( n ),              
    as.integer( K ),                   
    as.integer( m )                    # max leaves per node
  )
 
  #print( r.result[[2]] ) 
  #return( r.result[[2]] ) 

  return( (1:n)[ r.result[[2]] > .5 ] )
  
}



library(BalancedSampling) 

N <- 100000
n <- 1000
x <- cbind( runif(N), runif(N)) 


Cprog <- proc.time()
sampled2 <- lpm2( rep(n/N,N),x  )
print("lpm2 running time")
print(proc.time() - Cprog) 

Cprog <- proc.time()
sampled3 <- lpm3( rep(n/N,N),x  )
print("lpm3 running time")
print(proc.time() - Cprog) 

print(length(sampled3))

