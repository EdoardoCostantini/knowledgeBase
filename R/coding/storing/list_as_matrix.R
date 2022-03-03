# How smartly store in a matrix instead of list
# Subsetting guide 1

# Case: you generate many different datasets from a multivariate distribution
#       you wnat to generate the data, and then compute the covariance matrix
#       to each.

  n <- 15 # sample size
  p <- 5 # dimensionality
  N <- 50 # number of datasets
  
  Sigma_true <- diag(p)
  
  dataStorage <- matrix(rep(NA, n*p*N), ncol = p)

# Store data

  for(i in 1:N){
    
    row_indx <- ((i-1)*n + 1):( (i)*n )
    
    dataStorage[row_indx, ] <- mvtnorm::rmvnorm(n, rep(0, p), Sigma_true)
    
  }

# Compute covariance

  for(i in 1:N){
  
    row_indx <- ((i-1)*n + 1):( (i)*n )
    
    print(cov(dataStorage[row_indx, ]))
    
  }