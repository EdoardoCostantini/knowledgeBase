# Project:   knowledgeBase
# Objective: explore the ridge principle and ridge regression nuances
# Author:    Edoardo Costantini
# Notion:    https://lavish-hollyhock-981.notion.site/Ridge-regression-8134d8babda5413ab182df645c6196a8
# Other ref:
# Created:   2022-02-28
# Modified:  2022-02-28

# Centering X: effect on interpretation of OLS intercept -----------------------

# Take the mtcars data
y <- mtcars[, "mpg"]
X <- mtcars[, -1]

# Create a version of X that is centered
X_center <- scale(X, center = TRUE, scale = FALSE)

# Fit an regular linear model
lm_ols <- lm(y ~ X_center)

# Check that b0 is equal to the mean of y
coef(lm_ols)["(Intercept)"] - mean(y)

# Fitting ridge regression: manual and glmnet use ------------------------------

# Scale the data (standardize)
X_scale <- scale(X, center = TRUE, scale = TRUE)

# Compute the cross-product matrix of the data
XtX <- t(X_scale) %*% X_scale

# Define the identify matrix
I <- diag(ncol(X_scale))

# Define a lambda value
lambda <- .1

# Estimate the regression coefficients with the ridge penalty
betas <- solve(XtX + lambda * I) %*% t(X_scale) %*% y
