# Project:   knowledgeBase
# Objective: Explore a simple case of EM algorithm with missing values
# Author:    Edoardo Costantini
# Created:   2021-11-17
# Modified:  2022-03-01
# Notion:
# Other ref:

rm(list = ls())

# Packages ----------------------------------------------------------------

  library(mvtnorm)           # for likelihood functions
  library(ISR3)              # for SWP functions
  library(fastmatrix)        # alternative sweep
  library(mice)              # for NA pattern assesment
  library(rbenchmark)        # to benchmark code
  source("emSchafer.R") # main sweep function
  source("emFast.R")    # main sweep function

# Likelihood summary ------------------------------------------------------

# Review a few concepts you might need
  # Given some dataset (e.g. Little Rubin 2002 example 7.7, p. 152)
    Y_full <- matrix(c(7, 1, 11, 11, 7, 11, 3, 1, 2, 21, 1, 11, 10,
                       26, 29, 56, 31, 52, 55, 71 ,31, 54, 47,40,66,68,
                       6, 15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 9, 8,
                       60,52, 20, 47, 33, 22,6,44,22,26,34,12,12,
                       78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,93.1,
                       115.9,83.8,113.3,109.4), ncol = 5)
    n <- nrow(Y_full)
    mu <- colMeans(Y_full) # sample, but imagine it's population
    Sigma <- cov(Y_full)   # sample, but imagine it's population
    y_i_1ton <- mvtnorm::rmvnorm(n, mu, Sigma) # y_i ~ MN(mu, Sigma)

  # Single row likelihood
    lkl_rowi <- det(2*pi*Sigma)^(-1/2) * exp((-1/2)*t(Y_full[2, ] - mu) %*% solve(Sigma) %*% (Y_full[2, ] - mu)) # for i = 2
    mvtnorm::dmvnorm(Y_full[2,], mu, Sigma) # same
    log(lkl_rowi)

  # Complete-data likelihood
    log(det(Sigma)^(-n/2) * exp((-1/2)* sum(diag((t(t(Y_full)-mu) %*% solve(Sigma)) %*% (t(Y_full)-mu)))))
    # disregarding the proportionality constant
    lkl_comp <- sum(dmvnorm(Y_full, mu, Sigma))
    log(lkl_comp)

  # Maximum-likelihood estiamtes
    # Define sufficient statistics T1 and T2
    T1 <- colSums(Y_full)
    T2 <- t(Y_full) %*% Y_full # cross-product matrix

    # Re-write complete-data loglikelihood based on sufficient statistics
    - n/2*log(det(Sigma)) - n/2*t(mu)%*%solve(Sigma)%*%mu + t(mu)%*%solve(Sigma)%*%T1 - 1/2*sum(diag(solve(Sigma)%*%T2))

    # Estimates
    y_bar <- n^(-1) * T1
    S <- n^(-1) * (t(Y_full)-y_bar) %*% t(t(Y_full)-y_bar)

# Datasets ----------------------------------------------------------------

# Multivariate Monotone example 7.7 in Little Rubin 2002 (p.153) ####
  dat_mm <- matrix(data = c(7,1, 11, 11, 7, 11, 3, 1, 2, NA, NA, NA, NA,
                            26,29, 56, 31, 52, 55, 71 ,31, 54, NA,NA,NA,NA,
                            6,15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 9, 8,
                            60,52, 20, 47, 33, 22,NA,NA,NA,NA,NA,NA,NA,
                            78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,
                            93.1,115.9,83.8,113.3,109.4),
                   ncol = 5,
                   dimnames = list(NULL, c("X1", "X2", "X3", "X4", "X5"))
  )

  # NAs pattern
  md.pattern(dat_mm)

# Multivariate General pattern example 5.3.7 in Schafer 1997 (p.230) ####
  dat_mg <- matrix(data = c(16,20,16,20,-6,-4, 12,24,12,-6,4,-8,
                            8,8,26,-4,4,8, 20,8,NA,NA,20,-4,
                            8,4,-8,NA,22,-8, 10,20,28,-20,-4,-4,
                            4,28,24,12,8,18, -8,20,24,-3,8,-24,
                            NA,20,24,8,12,NA),
                   ncol = 6,
                   byrow = TRUE,
                   dimnames = list(NULL, c("15P","15L","15H", "90P", "90L", "90H")))

  # NAs pattern
  md.pattern(dat_mg)

# Prepare data ------------------------------------------------------------

  # EM Estimate of descriptives for data with missing values
  Y <- dat_mg # chose dataset here
  n <- nrow(Y)
  p <- ncol(Y)

  # Starting value
  theta0 <- matrix(rep(NA, (p+1)^2 ), ncol = (p+1),
                  dimnames = list(c("int", colnames(Y)),
                                  c("int", colnames(Y))
                  ))
  theta0[, 1]   <- c(-1, colMeans(Y, na.rm = TRUE)) # T1 CC
  theta0[1, ]   <- c(-1, colMeans(Y, na.rm = TRUE))
  theta0[-1,-1] <- cov(Y, use = "pairwise.complete.obs") * (n - 1)/n # T2 CC

  # Run function
  theta_hat <- emSchafer(Y = Y, iters = 500, theta0 = theta0)
  theta_hat_f <- emFast(Y = Y, iters = 500, theta0 = theta0)

  # Assess results
  theta0
  round(theta_hat[1, ], 3)             # means
  round(sqrt(diag(theta_hat)[-1]), 3)  # standard deviations
  round(theta_hat[-1, -1], 3)          # covariance matrix
  round(cov2cor(theta_hat[-1, -1]), 3) # correlation matrix

  # Compare speed
  round(theta_hat - theta_hat, 3)      # same results
  benchmark(
    "Scahfer" = { emSchafer(Y = Y, iters = 500, theta0 = theta0) },
    "Better R code" = { emFast(Y = Y, iters = 500, theta0 = theta0) },
    replications = 10,
    columns = c("test", "replications", "elapsed",
                "relative", "user.self", "sys.self")
  )

# RNG data ---------------------------------------------------------------------

  set.seed(20220327)
  n <- 1e4
  p <- 3
  Sig <- matrix(rep(1.5, p*p), ncol = p)
    diag(Sig) <- 3
  Y <- MASS::mvrnorm(n, rep(10, p), Sig)
  colnames(Y) <- paste0("y", 1:p)

# Estimate covmat on the original data

  cov_og <- cov(Y) * (n-1) / n # ML version
  muv_og <- colMeans(Y)

# Impose missingness

  amputeY <- mice::ampute(Y,
                          prop = .5,
                          patterns = matrix(c(1, 1, 0,
                                              1, 0, 1,
                                              1, 0, 0), ncol = 3, byrow = TRUE),
                          mech = "MAR",
                          type = "RIGHT"
  )

  Ymiss <- as.matrix(amputeY$amp)

# Estimate covmat on the listwise deletion data

  cov_cc <- cov(Ymiss, use = "complete.obs") * (n-1) / n # ML version
  muv_cc <- colMeans(Ymiss, na.rm = TRUE)

# Estimate covmat with EM
  # Give a bad starting value
  theta0 <- matrix(rep(NA, (p+1)^2 ), ncol = (p+1),
                   dimnames = list(c("int", colnames(Y)),
                                   c("int", colnames(Y))
                   ))
  theta0[, 1]   <- c(-1, rep(0, p)) # T1 CC
  theta0[1, ]   <- c(-1, rep(0, p))
  theta0[-1,-1] <- diag(p)

  # Run function
  theta_hat <- emSchafer(Y = Ymiss, iters = 500, theta0 = theta0)
  theta_hat_f <- emFast(Y = Ymiss, iters = 500, theta0 = theta0)

  # Extract results
  cov_EM <- theta_hat_f[-1, -1]
  muv_EM <- theta_hat_f[-1, 1]

# Compare results

  # function to simplify comparison between covmats
  vectomat <- function (comat){
          names_vec <- unlist(lapply(rownames(comat), function (i){
            paste0(i, "-", rownames(comat))
          }))
          comat_vec <- as.vector(comat)
          names(comat_vec) <- names_vec
          # Only unique
          return(comat_vec[!duplicated(comat_vec)])
        }

  # compare covmats
  cbind(OG = vectomat(cov_og),
        cc = vectomat(cov_cc),
        EM = vectomat(cov_EM),
        diff_cc = round(vectomat(cov_og) - vectomat(cov_cc), 3),
        diff_EM = round(vectomat(cov_og) - vectomat(cov_EM), 3)
  )


# Weighted version -------------------------------------------------------------

# use example dataset from example of stats::cov.wt()
Y <- Y_full <- cbind(Y1 = 1:10, Y2 = c(1:3, 8:5, 8:10))

# Add a missing value
Y[10, 1] <- NA
Y[1, 2] <- NA

p <- ncol(Y)
n <- nrow(Y)

# Starting value
theta0 <- matrix(rep(NA, (p+1)^2 ), ncol = (p+1),
                dimnames = list(c("int", colnames(Y)),
                                c("int", colnames(Y))
                ))
theta0[, 1]   <- c(-1, colMeans(Y, na.rm = TRUE)) # T1 CC
theta0[1, ]   <- c(-1, colMeans(Y, na.rm = TRUE))
theta0[-1,-1] <- cov(Y, use = "pairwise.complete.obs") * (n - 1)/n # T2 CC

# Run function
theta_hat_f <- emFast(Y = Y, iters = 1e3, theta0 = theta0)

# Target
cov.wt(Y, wt = w1, method = "ML")$cov # i.e. method = "unbiased"

# Define weights
w1 <- c(0, 0, 0, 1, 1, 1, 1, 1, 0, 0)
