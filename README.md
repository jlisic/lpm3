This is a "faster" implementation of the Local Pivotal Method for balanced sampling. 


  Current timings relative to BalancedSampling's lpm2 procedure on an i5-4258U CPU @ 2.40GHz running R 3.2.1 on OS X 10.10.4:

  N = 100,000
  n = 1,000
  x = cbind(1:N, 1:N)

  lpm2 (Balanced Sampling):
  run time = 43.708 elapsed seconds

  lpm3:
  run time = 0.183 elapsed seconds (Note, this includes the time to generate x, but does not include R function call and copying overhead). 

