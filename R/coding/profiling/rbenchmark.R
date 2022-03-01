# Project:   knowledgeBase
# Objective: benchmark profiling example
# Author:    Edoardo Costantini
# Created:   2021-01-29
# Modified:  2022-02-28
# Notion:
# Other ref:

rm(list = ls())

# Set up ------------------------------------------------------------------

  # Load the rbenchmark pacakge
  library(rbenchmark)

  # Define parameter for the MVN distribution
  nobs <- 4 # 1e3 # number of observations original datasets
  nvar <- 3 # 1e2 # number of columns original datasets
  npcs <- 2 # 5   # number of pcs extracted
  reps <- 3 # 1e3 # matrix multiplications

  # Define items and PC loadings
  X <- MASS::mvrnorm(nobs, rep(0, nvar), diag(nvar))
  P <- MASS::mvrnorm(nvar, rep(0, npcs), diag(npcs))

  # Create storing objects
  datasets <- lapply(1:reps, function(r) X)
  loadings <- lapply(1:reps, function(r) P)
  pairedList <- lapply(1:reps, function(x) list(X, P))

  # Define the function we want to benchmark
  matProd <- function(X, P){X %*% P}

# Benchmark --------------------------------------------------------------------
# Name all the different versions we are benchmarking within the statement

  benchmark(
    # use of mapply
    "mapply" = {
      out_mapply <- mapply(matProd, datasets, loadings, SIMPLIFY = FALSE)
    },
    # use of map
    "map" = {
      out_map <- Map(matProd, datasets, loadings)
    },
    # use of standard loop
    "loop" = {
      out_loop <- vector("list", reps)
      for(r in 1:reps){
        out_loop[[r]] <- datasets[[r]] %*% loadings[[r]]
      }
    },
    "lappy" = {
      lappyOut <- lapply(pairedList, function(x) Reduce("%*%", x))
    },
    # How many times we want to run the comparison
    replications = 1e3,
    # What are we interested in
    columns = c("test", "replications", "elapsed",
                "relative", "user.self", "sys.self")
  )
