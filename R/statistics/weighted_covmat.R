# Project:   knowledgeBase
# Objective: see how the weigthed covairance matrix work
# Author:    Edoardo Costantini
# Created:   2022-02-28
# Modified:  2022-03-02
# Notion:
# Other ref: https://en.wikipedia.org/wiki/Sample_mean_and_covariance#Weighted_samples

# Set up -----------------------------------------------------------------------

  # Get the dataset used in the example of stats::cov.wt()
  xy <- cbind(x = 1:10, y = c(1:3, 8:5, 8:10))
  w_in <- c(0,0,0,1,1,1,1,1,0,0)
    # weights must be non-negative and not all zero

  # Get the weighted estimate to define the target
  cov.wt(xy, wt = w_in) # i.e. method = "unbiased"

# cov.wt manually --------------------------------------------------------------

  # Assign values to the function arguments
  x      <- xy
  wt     <- w_in
  method <- "ML"

  # Assign values to some of the internal objects
  n <- nrow(x)
  p <- ncol(x)

  # Normalise weights (to sum to 1)
  wt <- wt / sum(wt)

  # Center on weighted mean if required
  center <- colSums(wt * x)

  # Center X on the weigthed mean
  x_cent <- sweep(x, 2, center, check.margin = FALSE)

    # Note that the sweep here is the same as subtracting the "center" to each value
    all.equal(
      sweep(x, 2, center, check.margin = FALSE),
      t(apply(x, 1, function (i) i - center))
    )

  # Weight data
  x_weighted <- sqrt(wt) * x_cent

  # compute cov based on standard cross-product (ML method)
  # The following are equivalent:
  cov.wt(xy, wt = w_in, method = "ML", center = TRUE)$cov
  crossprod(x_weighted)
  t(sqrt(wt) * x_cent) %*% (sqrt(wt)*x_cent)
  t(wt * x_cent) %*% (x_cent)
  t(w_in * x_cent) %*% (x_cent) / sum(w_in)

  # compute cov based on "corrected" cross-product (unbiased)
  # The following are equivalent:
  cov.wt(xy, wt = w_in, method = "unbiased", center = TRUE)$cov
  crossprod(x_weighted)/(1 - sum(wt^2))
  1 / (1 - sum(wt^2)) * t(wt * x_cent) %*% (x_cent)
  t(sqrt(wt) * x_cent) %*% (sqrt(wt)*x_cent) / (1 - sum(wt^2))
  t(wt * x_cent) %*% (x_cent) / (1 - sum(wt^2))
  t(w_in * x_cent) %*% (x_cent) / (1 - sum(w_in^2))

  # Tobs: the matrix of sufficient stats can be obtained using the unnormalized
  #       weights
  t(w_in * x_cent) %*% (x_cent)
  crossprod(sqrt(w_in) * x_cent)

# Long (loop) version ----------------------------------------------------------

  x      <- as.matrix(xy)
  wt     <- w_in
  n      <- nrow(x)
  s      <- sum(wt)
  wt     <- wt/s # normalized wieghts
  center <- colSums(wt * x)

  # Center data
  x_cent <- base::sweep(x, 2, center, check.margin = FALSE)

  # Compute the matrix of observed sufficient statistics (wegithed)
  Tobs <- matrix(0, ncol(x), ncol(x))
  for(i in 1:nrow(x)){
    Tobs <- Tobs + w_in[i] * (x_cent[i, ]) %*% t(x_cent[i, ])
  }

  # Compute the weighted sample size
  n_wt <- sum(w_in)

  # Compute the covariance matrix
  (covmat <- (Tobs / n_wt))
  cov.wt(xy, wt = wt, method = "ML")$cov