# Project:   knowledgeBase
# Objective: see how the weigthed covairance matrix work
# Author:    Edoardo Costantini
# Notion:
# Other ref:
# Created:   2022-02-28
# Modified:  2022-03-01

# Set up -----------------------------------------------------------------------

# Get the dataset used in the example of stats::cov.wt()
xy <- cbind(x = 1:10, y = c(1:3, 8:5, 8:10))
w1 <- c(0,0,0,1,1,1,1,1,0,0)
  # weights must be non-negative and not all zero

# Get the weighted estimate to define the target
cov.wt(xy, wt = w1) # i.e. method = "unbiased"

# Replicate manually -----------------------------------------------------------

# Assign values to the function arguments
x = xy
wt = w1
cor = FALSE
center = TRUE
method = "ML"

# Assign values to some of the internal objects
n <- nrow(x)

# Normalise weights (to sum to 1)
wt <- wt / sum(wt)

# Center on weighted mean if required
center <- colSums(wt * x)

# Transform X based on the weigthing
x <- sqrt(wt) * sweep(x, 2, center, check.margin = FALSE)

  # Note that the sweep here is the same as subtracting the "center" to each value
  all.equal(
    sweep(x, 2, center, check.margin = FALSE),
    t(apply(x, 1, function (i) i - center))
  )

# Compute the covariance matrix

# based on "corrected" cross-product (unbiased)
cov <- crossprod(x)/(1 - sum(wt^2))

# based on standard cross-product (ML method)
cov <- crossprod(x)

# Store the result
y <- list(cov = cov, center = center, n.obs = n)

# Output
y$cov
cov.wt(xy, wt = w1, method = "ML", center = FALSE)$cov