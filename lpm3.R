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
 
  n <- NROW(x)
  perm <- sample(0:(n-1),size=n) 
  p <- length(x) / n

  Cprog <- proc.time()

  # send our data to the C program
  r.result <- .C("R_lpm3",
    as.double( t(x) ),                 # 1 data we query
    as.double( prob ),                 # 2 probability
    as.integer( n ),              
    as.integer( p )            
  )
  
  print("C running time")
  print(proc.time() - Cprog) 

  return(r.result)
  
}




sampled <- lpm3( rep(1000/20000,20000), 1:20000 )




