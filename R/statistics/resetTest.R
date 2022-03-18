### Title:    Optional material: RESET test in more depth
### Author:   Kyle M. Lang, L.V.D.E. Vogelsmeier, Edo
### Created:  2022-03-18
### Modified: 2022-03-18
### Sources:  https://www.youtube.com/watch?v=DC7dAI4D-EI
### Tag:      MLR; Assumptions; non-linear    

# Set up -----------------------------------------------------------------------

# Load packages
  library(lmtest)

# Generate random predictors
  set.seed(1234)
  n <- 1e3
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  
# Generate a y with only linear effects
  y_lin <- 1 + x1 + x2 + x3 + rnorm(n)

# Fit linear model
  lm2 <- lm(y_quad ~ x1 + x2 + x3)

# Check RESET test
  resettest(lm2, power = 2, type = "fitted")    # H0 not rejected!
  resettest(lm2, power = 2, type = "regressor") # H0 not rejected!

# Generate a y with a quadratic term
  y_quad <- 1 + x1 + x2 + x3 + x3^2 + rnorm(n)
  
# Fit linear model
  lm2 <- lm(y_quad ~ x1 + x2 + x3)

# Perform RESET test
  resettest(lm2, power = 2, type = "fitted") # squared regressor is good
  
# Replicate the RESET test yourself --------------------------------------------

  # Fit linear model with no squared terms
  out1 <- lm(y_quad ~ x1 + x2 + x3)
  summary(out1)

  # Fit a lienar model predicting the dv with the squared fitted values from the original model
  out2 <- lm(y_quad ~ x1 + x2 + x3 + I(fitted(out1)^2))
  summary(out2)

  # Perform RESET w/ the R function
  RESET_fun <- resettest(out1,
                         power = 2,
                         type = "fitted") # using the squared regressors

  # Perform RESET test yourself
  RESET_asR2 <- anova(out1, out2) # It's a simple change in R-squared test

  # Same F statistic
  RESET_fun$statistic - RESET_asR2$F[2]
  
  # Same degrees of freedom
  RESET_regs$parameter - RESET_regs$parameter
