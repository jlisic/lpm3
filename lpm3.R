if ( !( 'loaded' %in% ls() ) )  {
  print('lpm3.so loaded')
  dyn.load('~/src/lpm3/lpm3.so')
  loaded <- T
}

  


############################################################
# x data set
# prob probability
############################################################
lpm3 <- function(
  prob,
  x                      # data to mode seek on
  ) {

 stop("not working yet")

  n <- NROW(x)
  p <- length(x) / n


  # send our data to the C program
  r.result <- .C("R_lpm3",
    as.double( t(x) ),                 # 1 data we query
    as.double( prob ),                 # 2 probability
    as.integer( n ),              
    as.integer( p )            
  )
  

  return(r.result)
  
}



library(BalancedSampling) 

N <- 100000
n <- 1000
x <- cbind( 1:N, 1:N) 


Cprog <- proc.time()
sampled <- lpm2( rep(n/N,N),x  )
print("lpm2 running time")
print(proc.time() - Cprog) 




