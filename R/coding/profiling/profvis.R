# Project:   knowledgeBase
# Objective: profvis profiling example
# Author:    Edoardo Costantini
# Notion:
# Other ref:
# Created:   2021-01-29
# Modified:  2022-02-28

rm(list = ls())

# Define parameter for the MVN distribution
n <- 1e3
p <- 1e3
rho <- .8
Sigma <- matrix(rho, nrow = p, ncol = p)
diag(Sigma) <- 1
mu <- rep(0, p)

# Run the two sampling functions you want to profile within a profvis statement
profvis::profvis({
  mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
})
