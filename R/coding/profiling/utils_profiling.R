# Project:   knowledgeBase
# Objective: utils profiling example
# Author:    Edoardo Costantini
# Notion:
# Other ref:
# Created:   2022-02-28
# Modified:  2022-02-28

rm(list = ls())

# Profile withing session ------------------------------------------------------

  # Open profiling
  Rprof()

  # Define parameter for the MVN distribution

  n <- 1e3
  p <- 1e3
  rho <- .8
  Sigma <- matrix(rho, nrow = p, ncol = p)
  diag(Sigma) <- 1
  mu <- rep(0, p)

  # Run the two sampling functions you want to profile

  mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)

  # Close profiling

  Rprof(NULL)
  out <- summaryRprof()
  out$by.total
  out$by.self

# Profile with costum filename -------------------------------------------------

  # Set up profiling
  Rprof(filename = "./utils_profiling.out")

  # Run the two sampling functions you want to profile
  mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)

  # Close profiling
  Rprof(NULL)
  out <- summaryRprof(filename = "./utils_profiling.out")

# Extract profile info on target functions -------------------------------------

  # Define target funcitons
  target <- c("mvtnorm::rmvnorm", "MASS::mvrnorm")
  rows <- paste0("\"", target, "\"")

  # Open profiling
  Rprof()

  # Run the two (target) sampling functions you want to profile
  mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)

  # Close profiling
  Rprof(NULL)

  # Extract what is of interest
  out <- summaryRprof()
  summaryRprof()$by.total[rows, "total.time"]

  # Profile within RStudio
  profvis::profvis({
    mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
    MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  })
