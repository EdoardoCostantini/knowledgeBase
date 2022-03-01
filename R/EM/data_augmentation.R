# Project:   knowledgeBase
# Objective: showing differences in how EM and DA missing data handling works
# Author:    Edoardo Costantini
# Created:   2021-11-17
# Modified:  2021-12-02
# Notion:
# Other ref:

# Packages --------------------------------------------------------
  library(mvtnorm) # for likelihood functions
  library(ISR3)    # for SWP functions
  library(mice)

# Functions Needed --------------------------------------------------------
  # Different functions that you need in the code to perform different things
  # to create the correct theta matrix
  # fill_cova <- function(dt, cova){
  #   n <- nrow(dt)
  #   for (j in 1:nrow(cova)) {
  #     for (l in 1:ncol(cova)) {
  #       if(j==l) cova[j,l] <- sum((dt[,j]-mean(dt[,j]))**2)/n
  #       if(j!=l) {
  #         cova[j,l] <- sum((dt[,j]-mean(dt[,j]))*(dt[,l]-mean(dt[,l])))/n
  #         cova[l,j] <- sum((dt[,j]-mean(dt[,j]))*(dt[,l]-mean(dt[,l])))/n
  #       }
  #     }
  #   }
  #   return(cova)
  # }
  # ad hoc variance covariance functions (biased versions)
  var_nden <- function(x){round(sum( ((x-mean(x))**2) )/length(x), 3)}
  cov_nden <- function(x,y){round(sum( ((x-mean(x))*(y-mean(y))) )/length(x), 3)}

  # to add matrices stored as elements of a list
  add_matrices <- function(listMat){
    summedMat <- listMat[[1]]
    for (m in 2:length(listMat)) {
      summedMat <- summedMat+listMat[[m]]
    }
    return(summedMat)
  }

# Likelihood summary ------------------------------------------------------
# Review a few concepts you might need
  # Given some dataset (e.g. Little Rubin 2002 ecxample 7.7 full data)
    Y_full <- matrix(c(7,1, 11, 11, 7, 11, 3, 1, 2, 21, 1, 11, 10,
                       26,29, 56, 31, 52, 55, 71 ,31, 54, 47,40,66,68,
                       6,15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 9, 8,
                       60,52, 20, 47, 33, 22,6,44,22,26,34,12,12,
                       78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,93.1,
                       115.9,83.8,113.3,109.4), ncol = 5)
    n <- nrow(Y_full)
    mu <- colMeans(Y_full) # sample, but imagine it's population
    Sigma <- cov(Y_full)   # sample, but imagine it's population
    y_i_1ton <- mvtnorm::rmvnorm(n, mu, Sigma) # y_i ~ MN(mu, Sigma)

  # Single row likelihood
    lkl_rowi <- det(2*pi*Sigma)**(-1/2) * exp((-1/2)*t(Y_full[2, ] - mu) %*% solve(Sigma) %*% (Y_full[2, ] - mu)) # for i = 2
      dmvnorm(Y_full[2,], mu, Sigma) # same
    log(lkl_rowi)

  # Complete-data likelihood
    log(det(Sigma)**(-n/2) * exp((-1/2)* sum(diag((t(t(Y_full)-mu) %*% solve(Sigma)) %*% (t(Y_full)-mu)))))
    # disregarding the proportionality constant
    lkl_comp <- sum(dmvnorm(Y_full, mu, Sigma))
    log(lkl_comp)

  # Maximum-likelihood estiamtes
    # Define sufficient statistics T1 and T2
    T1 <- as.vector(t(Y_full) %*% rep(1, nrow(Y_full))) # vector of column sums (sum of values of each variable)
    colSums(Y_full)

    T2 <- t(Y_full) %*% Y_full # matrix of columnwise sums of squares (diagonals) and crossproducts (offdiagonals)
      Y_full[,1] %*% Y_full[,1]  # repeat for all columns -> T2

    # Re-write complete-data loglikelihood
    - n/2*log(det(Sigma)) - n/2*t(mu)%*%solve(Sigma)%*%mu + t(mu)%*%solve(Sigma)%*%T1 - 1/2*sum(diag(solve(Sigma)%*%T2))

    # Estimates
    y_bar <- n**(-1) * colSums(Y_full)
    S <- n**(-1) * (t(Y_full)-y_bar) %*% t(t(Y_full)-y_bar)

# . -----------------------------------------------------------------------
# Datasets ----------------------------------------------------------------
# List of datasets with different caracteristics so that you can explore
# different things.
# > Univariate Example 5.3.6 (Schafer 1997, p.226) ####
  Y <- matrix(c(270,236,210,142,280,272,160,220,226,242,186,266,206,318,294,282,234,224,276,282,360,310,280,278,288,288,244,236,
                218,234,214,116,200,276,146,182,238,288,190,236,244,258,240,294,220,200,220,186,352,202,218,248,278,248,270,242,
                156,NA,242,NA,NA,256,142,216,248,NA,168,236,NA,200,264,NA,264,NA,188,182,294,214,NA,198,NA,256,280,204), ncol = 3)
  colnames(Y) <- c("Y1","Y2","Y3")
  p <- ncol(Y)

  # Starting values
  theta <- matrix(rep(NA, (p+1)*(p+1) ), ncol = (p+1))
  theta[,1] <- c(1, 200, 200, 200) # variable means in the first column
  theta[1,] <- c(1, 200, 200, 200) # (and first row)
  theta[-1,-1] <- 50**2*diag(p)
  theta0 <- theta

# > Multivariate monotone Little Rubin 2002 example 7.7 ####
  Y <- matrix(c(7,1, 11, 11, 7, 11, 3, 1, 2, NA, NA, NA, NA,
                 26,29, 56, 31, 52, 55, 71 ,31, 54, NA,NA,NA,NA,
                 6,15, 8, 8, 6, 9, 17, 22, 18, 4, 23, 9, 8,
                 60,52, 20, 47, 33, 22,NA,NA,NA,NA,NA,NA,NA,
                 78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,93.1,115.9,83.8,113.3,109.4), ncol = 5)
  colnames(Y) <- c("X1", "X2", "X3", "X4", "X5")
  p <- ncol(Y)
  n <- nrow(Y)
  # Y <- Y[, c(3,5,1,2,4)] # monotonly oredered data
  # Starting values
  theta <- matrix(rep(NA, (p+1)*(p+1) ), ncol = (p+1))
  theta[,1] <- c(-1, colMeans(Y, na.rm = TRUE)) # variable means in the first column
  theta[1,] <- c(-1, colMeans(Y, na.rm = TRUE))
  theta[-1,-1] <- cov(Y, use = "complete.obs")
  theta0 <- theta

# > Multivariate example 5.3.7 (Schafer 1997, p.230) ####
  Y <- matrix(c(16,20,16,20,-6,-4,
                12,24,12,-6,4,-8,
                8,8,26,-4,4,8,
                20,8,NA,NA,20,-4,
                8,4,-8,NA,22,-8,
                10,20,28,-20,-4,-4,
                4,28,24,12,8,18,
                -8,20,24,-3,8,-24,
                NA,20,24,8,12,NA), ncol = 6, byrow = TRUE)
  colnames(Y) <- c("15P","15L","15H", "90P", "90L", "90H")
  p <- ncol(Y)

  # Starting values
  theta <- matrix(rep(NA, (p+1)*(p+1) ), ncol = (p+1))
  theta[,1] <- c(-1, colMeans(Y, na.rm = TRUE)) # variable means in the first column
  theta[1,] <- c(-1, colMeans(Y, na.rm = TRUE))
  theta[-1,-1] <- cov(Y, use = "complete.obs")
  theta0 <- theta

# > Multivariate own dataset ####
  set.seed(20200128)
  p <- 10; n <- 100
  Sigma_true <- MCMCpack::riwish(p, diag(p))
  mu_true <- round(runif(p, 0, 100), 0)
  dt_comp <- mvtnorm::rmvnorm(n, mu_true, Sigma_true)
   colnames(dt_comp) <- paste0("Y", 1:ncol(dt_comp))

  missing_reason_1 <- dt_comp[,5] > mu_true[5]+sqrt(var(dt_comp[,5]))
  missing_reason_2 <- dt_comp[,5] < mu_true[5]-sqrt(var(dt_comp[,5]))
  missing_reason_4 <- dt_comp[,3] > mu_true[3]+sqrt(var(dt_comp[,3]))
  Y <- dt_comp
  Y[missing_reason_1, 1] <- NA
  Y[missing_reason_2, 2] <- NA
  Y[missing_reason_4, 4] <- NA

  mice::md.pattern(as.data.frame(Y))

  # Starting values
  theta <- matrix(rep(NA, (p+1)*(p+1) ), ncol = (p+1))
  theta[,1] <- c(1, colMeans(Y, na.rm = TRUE)) # variable means in the first column
  theta[1,] <- c(1, colMeans(Y, na.rm = TRUE))
  theta[-1,-1] <- cov(Y, use = "complete.obs")
  theta0 <- theta

# > Employee data Enders 2010 example 4.14 ####
  Y <- as.matrix(read.table("./data/employee.dat"))
  p <- ncol(Y)
  Y[Y==c(-99)] <- NA

  md.pattern(Y)

  # Starting values
  theta <- matrix(rep(NA, (p+1)*(p+1) ), ncol = (p+1))
  theta[,1] <- c(-1, colMeans(Y, na.rm = TRUE)) # variable means in the first column
  theta[1,] <- c(-1, colMeans(Y, na.rm = TRUE))
  theta[-1,-1] <- cov(Y, use = "complete.obs")
  theta0 <- theta

# . -----------------------------------------------------------------------
# Preparation of dataset --------------------------------------------------

  # Define Missing data patterns
  patts <- mice::md.pattern(Y)
  R <- patts[-nrow(patts),-ncol(patts)] # missing data patterns (1=variable observed in a specific pattern)
  R <- R[, colnames(Y)]
  S <- nrow(R) # number of missing data patterns
  O <- apply(R, 1, function(x) {colnames(R)[x==1]}) # subset of column lables corresponding to var observed for a given pattern
  M <- apply(R, 1, function(x) {colnames(R)[x==0]})

  # Define I matrices (which obs in which pattern)
  n <- nrow(Y)        # number of observations in Y
  tf_pat_type <- R==1 # pattern config saved as True and False
  I <- vector("list", S)
  for (s in 1:S) {
    # s <- 1
    index <- NULL
    for (i in 1:n) {
      # i <- 1
      if(all.equal(!is.na(Y)[i, ], tf_pat_type[s,]) == TRUE) {
        index <- c(index, i)
      }
    }
    I[[s]]<- index
  }

  # Define vector of var names
  vars <- colnames(Y)

  # Define sufficient statistics matrix (observed)
  Tobs_s <- vector("list", S)

  for (s in 1:S) {
    indx_obsi  <- I[[s]]
    indx_obsV <- (which(vars %in% O[[s]])) # column number of variables observed in pattern s
    indx_misV <- (which(vars %in% M[[s]])) # column number of variables missing in pattern s

    Tobs_s[[s]] <- matrix(rep(NA, (p+1)**2), ncol = (p+1))

    # Fill number of observations in this pattern
    Tobs_s[[s]][1, 1] <- length(I[[s]])

    # Fill sums of observed values and cross products for variables observed in this pattern
    Tobs_s[[s]][indx_obsV+1, indx_obsV+1] <- t(Y[indx_obsi, indx_obsV]) %*% (Y[indx_obsi, indx_obsV]) # sum of products

    if(length(indx_obsi)==1){
      Tobs_s[[s]][1, indx_obsV+1] <- Y[indx_obsi, indx_obsV] # the only value you have
      Tobs_s[[s]][indx_obsV+1, 1] <- Y[indx_obsi, indx_obsV]
    } else {
      Tobs_s[[s]][1, indx_obsV+1] <- colSums(Y[indx_obsi, indx_obsV]) # sum of observed values
      Tobs_s[[s]][indx_obsV+1, 1] <- colSums(Y[indx_obsi, indx_obsV])
    }

    # Fill 0 for variances and cross porducts with missing variables for this pattern
    Tobs_s[[s]][indx_misV+1, ] <- 0
    Tobs_s[[s]][, indx_misV+1] <- 0

    # Make it upper triangular (useful for understanding the EM code)
    # Tobs_s[[s]][lower.tri(Tobs_s[[s]])] <- NA
  }

# . -----------------------------------------------------------------------
# A1 - Algorithm w/ little use of sweeping and full explanation ----------------

  # Complete data estimates
  theta0[1, -1]
  mu_true[c(1, 2, 4)]

  # Reset algorithm
  theta <- theta0

for (iter in 1:20) {

  # One Iteration
  T <- Reduce("+", Tobs_s)

  #> E-step ####
  # Description: for an assumed value of theta, computes the complete-data sufficient statistics
  #              over P(Y_mis | Y_obs, theta). We can rewrite:
  #
  #                 P(Y_mis | Y_obs, theta) = product( P(y_i_mis | y_i_obs, theta) )
  for(s in 2:S){
    # S; s <- 3
    obs <- I[[s]]
    Y_obs_name <- O[[s]]
    Y_mis_name <- M[[s]]

    for(j in 1:p){
      # j=1
      if(R[s, j] == 1 & theta[j+1,j+1]>0) {theta <- SWP(theta, (j+1))}
      # Makes the current j-th variable a predictor in the multivariate regression of the following variables
      # on the preceeding variables.
      # e.g. when j = 1, the lower left part theta contains the regression coefficinets of Y2, Y3, Y4, Y5 ~ Y1 (or Y2 ~ Y1, Y3 ~ Y1, ..., Yp ~ Y1)
      #      when j = 3, the lower left part theta contains the regression coefficinets of Y4, Y5 ~ Y1 + Y2 + Y3
      #      when j = p, the lower left part theta contains the regression coefficinets of Yp ~ Y1 + Y2 + Y3 + ... + Y(p-1)
      if(R[s, j] == 0 & theta[j+1,j+1]<0) {theta <- RSWP(theta, c(j+1))} # book keeping: diagonal element is negative if and only if thta has been swept on position j
      # Basically, thanks to this line, a pattern that is fully observed does not influence the theta becuase this
      # line reverses any sweept that was based on that fully observed pattern.
    }

    for(i in 1:length( I[[s]] )){ # create expected value
      # i <- 1 # for the second individual
      iID <- I[[s]][i]
      cj <- NULL # empty vector to store y_ij_star (as many elements as there are missing variables in pattern s)

      print(paste("Individual:", iID, "| Missing data pattern:", s))

      for(j in 1:length( M[[s]] ) ){ # for each variable missing a value in this missing patter
        # length( M[[s]] ); j <-1 # for the first variable missing
        J <- which(vars == Y_mis_name[j])

        cj[j] <- theta[1,J+1] # a_0j, intercept for missing variable j of M(s)

        print(paste("1. For var w/ miss:", Y_mis_name[j]))

        print(paste("   1.a) Intercept:", round(cj[j],3) ))

        for(k in 1:length( O[[s]] )){
          # k <- 2
          K <- which(vars == Y_obs_name[k])
          b <- theta[K+1,J+1]          # regression coefficient for observed variables K in predicting missing variable J

          print(paste("   1.b) add:", round(b,3), "*", Y_obs_name[k]))

          cj[j] <- cj[j] + b * Y[iID,K] # predicted value for variable j (missing) for individual i (based on observed K)

        }
          print(paste("   Result: predicted value for individual", iID, ":", round(cj[j],3) ))
      }

      for(j in 1:length( M[[s]] ) ){ # update T matrix
        # length( M[[s]] ); j <- 1
        J <- which(vars == Y_mis_name[j])
        print(paste("2. For var w/ miss:", Y_mis_name[j]))

        print(paste("   add predicted value", round(cj[j],3), "to T1 at", J, "th spot"))

        T[1,J+1] <- T[1,J+1] + cj[j]
        # T[J+1,1] <- T[J+1,1] + cj[j]

        for(k in 1:length( O[[s]] )){
          # length( O[[s]] ); k <- 1
          K <- which(vars == Y_obs_name[k])

          # print(paste("   add cross product of observed", Y_obs_name[k], "(", Y[iID,K], ")",
          #             "and predicted", round(cj[j],3), "to T2"))

          print(paste("   add cross product of predicted", Y_mis_name[j], round(cj[j],3),
                      "and observed", Y_obs_name[k], Y[iID,K], "to T2 for theta", round(theta[K+1,J+1],3) ) )

          T[K+1,J+1] <- T[K+1,J+1] + cj[j]*Y[iID,K]
          # T[J+1,K+1] <- T[J+1,K+1] + cj[j]*Y[iID,K]
        }
        for(k in 1:length( M[[s]] )){
          # length( M[[s]] ); k <- 1
          K <- which(vars == Y_mis_name[k])

          if(K>=J){
            print(paste("   add cross product of predicted", Y_mis_name[j], round(cj[j],3),
                        "and predicted", Y_mis_name[k], round(cj[k],3), "to T2 for theta", round(theta[K+1,J+1],3)))
            T[K+1,J+1] <- T[K+1,J+1] + theta[K+1,J+1] + cj[k]*cj[j]
            # T[J+1,K+1] <- T[J+1,K+1] + theta[J+1,K+1] + cj[k]*cj[j]
            }
        }
      }
    print("--------------------------")
    print("----- NEW INDIVIDUAL -----")
    print("--------------------------")
    }
  }
  # Make T symmetric
    tT <- t(T)
    tT[upper.tri(T)] <- T[upper.tri(T)]
    T <- tT

  # > M-step ####
    theta <- SWP((n**(-1) * T), 1)

  # Print updated theta
    print(
      list(
        theta = theta,
        SDs = sqrt(diag(theta[-1,-1])),
        mu_hat = round(theta[-1,1],1)
        # sigma = sigma,
        # rho13 = theta[2,4]/(sigma[1,1]*sigma[3,3]), # rho13
        # rho23 = theta[3,4]/(sigma[2,2]*sigma[3,3]) # rho 23
      )
    )

}

# . -----------------------------------------------------------------------
# A2 - Algorithm w/ Sweeping shortcuts ------------------------------------

  # Reset algorithm
  theta <- theta0

  # One Iteration
  T <- add_matrices(Tobs_s)

  #> E-step ####
  for(s in 2:S){
    # S; s <- 2
    obs <- I[[s]]
    Y_obs_name <- O[[s]]
    Y_mis_name <- M[[s]]

    ####

    # REGULAR WAY
    #
    # for(j in 1:p){
    #   if(R[s, j] == 1 & thetaL[j+1,j+1]>0) {thetaL <- SWP(thetaL, (j+1))}
    #   # Makes the current j-th variable a predictor in the multivariate regression of the following variables
    #   # on the preceeding variables.
    #   # e.g. when j = 1, the lower left part theta contains the regression coefficinets of Y2, Y3, Y4, Y5 ~ Y1 (or Y2 ~ Y1, Y3 ~ Y1, ..., Yp ~ Y1)
    #   #      when j = 3, the lower left part theta contains the regression coefficinets of Y4, Y5 ~ Y1 + Y2 + Y3
    #   #      when j = p, the lower left part theta contains the regression coefficinets of Yp ~ Y1 + Y2 + Y3 + ... + Y(p-1)
    #   if(R[s, j] == 0 & thetaL[j+1,j+1]<0) {thetaL <- RSWP(thetaL, c(j+1))} # book keeping: diagonal element is negative if and only if thta has been swept on position j
    # }

    ####

    # SWEEP SHORTCUT (see also end of s loop for reset of theta)
    sweep_over <- which(vars %in% O[[s]])
    DVset <- which(vars %in% M[[s]])

    theta <- SWP(theta, c(sweep_over+1) )

    Bt_row_index <- c(DVset+1) # +1 accomodates the computer indexing (regular: 0:p -> computer: 0+1:p+1)
    Bt_col_index <- c(1, sweep_over+1)

    Bt <- theta[Bt_row_index, Bt_col_index] # lower left or upper right partition (reg coefs)
    C <- theta[(DVset+1), (DVset+1)]        # lower right partition (res variance)

    for(i in 1:length( I[[s]] )){ # create expected value
      # i <- 1
      ## Computing expecatitons ##

      iID <- I[[s]][i]
      cj <- NULL # empty vector to store y_ij_star (as many elements as there are missing variables in pattern s)

      for(j in 1:length( M[[s]] ) ){
        # length( M[[s]] ); j <- 1

        ####

        # REGULAR WAY

        # J <- which(vars == Y_mis_name[j])
        #
        # cj[j] <- theta[1,J+1] # a_0j, intercept for missing variable j of M(s)
        #
        # for(k in 1:length( O[[s]] )){
        #   # k <- 2
        #   K <- which(vars == Y_obs_name[k])
        #   b <- theta[K+1,J+1]        # regression coefficient for observed variables K in predicting missing variable J
        #   cj[j] = cj[j] + b * Y[iID,K] # dataset value for individual i, observed in variable K
        # }

        ####

        # SWEEP SHORTCUT
        if(is.matrix(Bt)==TRUE){ # If the pattern has more then 1 variables w/ missing values
          cj[j] <- Bt[j,] %*% as.vector(c(1, Y[iID, sweep_over])) # one regression for each variable w/ missing value
        }
        if(is.matrix(Bt)==FALSE){ # the pattern has more 1 variables w/ missing values (Bt is vector, not matrix in such a case)
          cj[j] <- Bt %*% as.vector(c(1, Y[iID, sweep_over])) # one regression for each variable w/ missing value
        }
      }

      ## updating T matrix ##

      # 1. adding expectation of y_ij (for mean sufficient statistics)

      T[1, DVset+1] <- T[1, DVset+1] + cj
      # T[DVset+1, 1] <- T[DVset+1, 1] + cj # expectation of yij is added to T1 (1st row and column of T)

      # 2. adding expectation of y_ij * y_ik (for covariance and variance sufficient statistics)

      for(j in 1:length( M[[s]] ) ){
        # j <- 1
        J <- which(vars == Y_mis_name[j])
        for(k in 1:length( O[[s]] )){
          # length( O[[s]] ); k <- 1
          K <- which(vars == Y_obs_name[k])
          T[K+1,J+1] <- T[K+1, J+1] + cj[j]*Y[iID,K] # adding expectation of y_ij * y_ik when
                                                     # k is observed for this iID
        }
        for(k in 1:length( M[[s]] )){
          # length( M[[s]] ); k <- 1
          K <- which(vars == Y_mis_name[k])
          if(K>=J){
            T[K+1,J+1] <- T[K+1,J+1] + theta[K+1,J+1] + cj[j]*cj[k] # adding expectation of y_ij * y_ik when
                                                                    # both j and k are unobserved for this iID
          }
        }
      }
    }
      theta <- RSWP(theta,  c(sweep_over+1) )
        # Book keeping: this corresponds to the reverse sweep in the RSWP loop performed in the
        # algorithm proposed by Schafer 1997. For one E step, the covariance matrix used to compute
        # individual contirbutions in each missing data pattern is the same! So after saving the
        # regression parameters needed for the computation of the expectations in one specifc
        # missing data pattern, I revert the theta to the one that every missing data pattern should
        # use at the current iteration. At the end of this iteration, the M step updates theta, and
        # all the missing data patterns in the following iteration will use the updated theta.
  }
  # Make T symmetric
    tT <- t(T)
    tT[upper.tri(T)] <- T[upper.tri(T)]
    T <- tT

  # > M-step ####
    theta <- SWP((n**(-1) * T), 1)

  # Print updated theta
    theta

# . -----------------------------------------------------------------------
# A3 - Algorithm functioning ----------------------------------------------
# Reset algorithm
theta <- theta0

# Iterate tot number of times
for (iter in 1:10) {

#> E-Step ####

  T <- add_matrices(Tobs_s) # start from observed data
  # add_matrices(Tobs_s) - T # good to inspect differences
  for (s in 2:S) {
    # S; s <- 3
    obs <- I[[s]]
    Y_obs_name <- O[[s]]
    Y_mis_name <- M[[s]]

    sweep_over <- which(vars %in% O[[s]])
    DVset <- which(vars %in% M[[s]])

    theta <- SWP(theta, c(sweep_over+1) )

    Bt_row_index <- c(DVset+1) # +1 accomodates the computer indexing (regular: 0:p -> computer: 0+1:p+1)
    Bt_col_index <- c(1, sweep_over+1)

    Bt <- theta[Bt_row_index, Bt_col_index] # lower left or upper right partition (reg coefs)
    C <- theta[(DVset+1), (DVset+1)]        # lower right partition (res variance)

    if(length(M[[s]])!=0) { # Perform only for patterns actually with missing variables

    for(i in 1:length( I[[s]] )){
      # length( I[[s]] ); i <- 1

      ## Computing expecatitons ##

      iID <- I[[s]][i]
      cj <- NULL # empty vector to store y_ij_star (as many elements as there are missing variables in pattern s)

      for(j in 1:length( M[[s]] ) ){
        # j <- 1
        if(is.matrix(Bt)==TRUE){
          cj[j] <- Bt[j,] %*% as.vector(c(1, Y[iID, sweep_over])) # one regression for each variable w/ missing value
        }
        if(is.matrix(Bt)==FALSE){
          cj[j] <- Bt %*% as.vector(c(1, Y[iID, sweep_over])) # one regression for each variable w/ missing value
        }
      }

      ## updating T matrix ##

      # 1. adding expectation of y_ij (for mean sufficient statistics)

      T[1, DVset+1] <- T[1, DVset+1] + cj # expectation of yij is added to T

      # 2. adding expectation of y_ij * y_ik (for covariance and variance sufficient statistics)

      for(j in 1:length( M[[s]] ) ){
        # length( M[[s]] ); j <- 1
        J <- which(vars == Y_mis_name[j])

        for(k in 1:length( O[[s]] )){
          # length( O[[s]] ); k <- 2
          K <- which(vars == Y_obs_name[k])
          T[K+1,J+1] <- T[K+1, J+1] + cj[j]*Y[iID,K] # adding expectation of y_ij * y_ik when
          # k is observed for this observation
        }

        for(k in 1:length( M[[s]] )){
          # length( M[[s]] ); k <- 1
          K <- which(vars == Y_mis_name[k])
          if(K>=J){
            T[K+1,J+1] <- T[K+1,J+1] + theta[K+1,J+1] + cj[k]*cj[j] # adding expectation of y_ij * y_ik when
            # both j and k are unobserved for this observation
          }
        }

      }

    }
    }

    # Reverse Sweep
    theta <- RSWP(theta,  c(sweep_over+1) )
      # End of pattern, revert the theta to the original one for this itarition so
      # so that every pattern starts from the same theta

  }
  # Make T symmetric
  tT <- t(T)
  tT[upper.tri(T)] <- T[upper.tri(T)]
  T <- tT

  # > M-step ####

  theta <- SWP((n**(-1) * T), 1)

  # Summaries iterations
  # sigma <- sqrt(theta[-1,-1])
  print(
    list(
      theta = theta,
      mu_hat = round(theta[-1,1],1)
      # sigma = sigma,
      # rho13 = theta[2,4]/(sigma[1,1]*sigma[3,3]), # rho13
      # rho23 = theta[3,4]/(sigma[2,2]*sigma[3,3]) # rho 23
    )
  )

}

# . -----------------------------------------------------------------------
# A4 - Data Augmentation Algorithm -----------------------------------
# Reset algorithm
theta <- theta0

# Iterations
  IT <- 1e3
  BURN <- 1e2
# Store posterior
  mu_GSigma <- matrix(rep(NA, (1e3+1e2)*p), ncol = p)
  Sigma <- vector("list", (1e3+1e2))

# Iterate tot number of times
for (iter in 1:(1e3+1e2)) {

  #> E-Step ####

  T <- add_matrices(Tobs_s) # start from observed data

  for (s in 2:S) {
    # S; s <- 2
    obs <- I[[s]]
    Y_obs_name <- O[[s]]
    Y_mis_name <- M[[s]]

    sweep_over <- which(vars %in% O[[s]])
    DVset <- which(vars %in% M[[s]])

    # Sweeping
    for(j in 1:p){
      if(R[s, j] == 1 & theta[j+1,j+1]>0) {theta <- SWP(theta, (j+1))}

      # Makes the current j-th variable a predictor in the multivariate regression of the following variables
      # on the preceeding variables.
      # e.g. when j = 1, the lower left part theta contains the regression coefficinets of Y2, Y3, Y4, Y5 ~ Y1 (or Y2 ~ Y1, Y3 ~ Y1, ..., Yp ~ Y1)
      #      when j = 3, the lower left part theta contains the regression coefficinets of Y4, Y5 ~ Y1 + Y2 + Y3
      #      when j = p, the lower left part theta contains the regression coefficinets of Yp ~ Y1 + Y2 + Y3 + ... + Y(p-1)
      if(R[s, j] == 0 & theta[j+1,j+1]<0) {theta <- RSWP(theta, c(j+1))} # book keeping: diagonal element is negative if and only if thta has been swept on position j
    }

    C <- chol(theta[DVset+1, DVset+1])

    if(length(M[[s]])!=0) { # Perform only for patterns actually with missing variables

      for(i in 1:length( I[[s]] )){ # create expected value
        # i <- 1
        ## Computing expecatitons ##

        iID <- I[[s]][i]
        # cj <- NULL # empty vector to store y_ij_star (as many elements as there are missing variables in pattern s)

        z <- NULL
        y_ij <- NULL
        for(j in 1:length( M[[s]] ) ){
          # length( M[[s]] ); j <- 3
          J <- which(vars == Y_mis_name[j])
          y_ij[j] <- theta[0+1, J+1]
          for(k in 1:length( O[[s]] )){
            # length( O[[s]] ); k <- 2
            K <- which(vars == Y_obs_name[k])
            y_ij[j] <- y_ij[j] + theta[K+1, J+1] * Y[iID,K]
          }
          z[j] <- rnorm(1,0,1)

          for(k in 1:length( M[[s]] )){
            # length( M[[s]] ); k <- 3
            K <- which(vars == Y_mis_name[k])
            if(K<=J){
              y_ij[j] <- y_ij[j] + C[k, j] * z[k]
            }
          }

        ## updating T matrix ##
        # used to to describe posteriors in P step

        # 1. adding expectation of y_ij (for mean sufficient statistics)
          T[0+1, J+1] = T[0+1, J+1] + y_ij[j]
          T[J+1, 0+1] = T[J+1, 0+1] + y_ij[j]

        # 2. adding expectation of y_ij * y_ik (for covariance and variance sufficient statistics)
          for(k in 1:length( O[[s]] )){
            # length( O[[s]] ); k <- 2
            K <- which(vars == Y_obs_name[k])
            T[K+1,J+1] <- T[K+1, J+1] + y_ij[j]*Y[iID,K] # adding expectation of y_ij * y_ik when
                                                         # k is observed for this observation
            T[J+1,K+1] <- T[J+1,K+1]  + y_ij[j]*Y[iID,K]
          }
          for(k in 1:length( M[[s]] )){
            # length( M[[s]] ); k <- 3
            K <- which(vars == Y_mis_name[k])
            if(K<=J){
              T[K+1,J+1] <- T[K+1, J+1] + y_ij[j]*y_ij[k] # Y[iID,K] # adding expectation of y_ij * y_ik when
                                                          # both j and k are NOT observed for this observation
              T[J+1,K+1] <- T[J+1,K+1]  + y_ij[j]*y_ij[k]
            }
          }
        # New variable j belonging to missing variables
        }
      # New individual i belonging to pattern s
      }
    }
  # New pattern s
  }

  # > M-step ####
  y_bar <- T[-1,1]/n

  # Hyperparameters updates
  # n <- nrow(Y)
  # m <- p+1 # m >= p
  # Lambda
  # tau <- 1
  # tau_prime <- tau + n
  # m_prime <- m + n
  # mu0_prime <-  (n/(t+n))*y_bar+(n/(t+n))*mu0
  #S <- n**(-1) * t(Y) %*% Y - y_bar %*% t(y_bar)
  nS <- T[-1,-1] # nS meaning n times S
  sqrt(diag(nS)/n)
  # Lambda_prime <- solve( solve(Lambda) + n * S + (tau*n / (tau+n)) * (y_bar-mu0)%*%t(y_bar-mu0) )
  #
  # Sigma <- MCMCpack::riwish(m, Delta)
  # mu_GSigma <- mvtnorm::rmvnorm(p, mu0, tau**(-1)*Delta)

  # Posterior draws using improper uninformative prior 5.18
  Sigma[[iter]] <- MCMCpack::riwish( n-1, nS/n)
  mu_GSigma[iter, ] <- mvtnorm::rmvnorm(1, y_bar, n**(-1)*Sigma[[iter]])

  theta[-1,1] <- mu_GSigma[iter, ]
  theta[1,-1] <- mu_GSigma[iter, ]
  theta[-1,-1] <- Sigma[[iter]]

# New iteration iter
}

  hist(mu_GSigma[-c(1:100),3])

  # Summaries iterations
  sigma <- sqrt(theta[-1,-1])
  print(
    list(
      theta = theta,
      mu_hat = round(theta[-1,1],1),
      sigma = sigma,
      rho13 = theta[2,4]/(sigma[1,1]*sigma[3,3]), # rho13
      rho23 = theta[3,4]/(sigma[2,2]*sigma[3,3]) # rho 23
    )
  )
