# Project:   knowledgeBase
# Objective: Explore how cell means coding effects computation of R-square
# Author:    Edoardo Costantini
# Created:   2022-03-03
# Modified:  2022-03-03
# Notion:
# Other ref:

# Set up -----------------------------------------------------------------------

  # Load the iris data
  data(iris)
  n <- nrow(iris)

# Default R-2 computation by summary() -----------------------------------------

  # Fit lm with dummy codes
  out_int <- lm(Petal.Length ~ Sepal.Length + Sepal.Width, data = iris)
  out_noi <- lm(Petal.Length ~ - 1 + Sepal.Length + Sepal.Width, data = iris)

  lm_fit <- out_noi
  r <- lm_fit$residuals
  f <- lm_fit$fitted.values
  rdf <- lm_fit$df.residual

  # Compute Mean squared error based on intercept presence
  if (attr(lm_fit$terms, "intercept")){
    # If the intercept was estimated, compute as:
    mss <- sum((f - mean(f))^2)
  } else {
    # If the intercept was not estimated, compute as:
    mss <- sum(f^2)
    # Note: this is correct only if the dependent variable has been centered on its mean
  }

  # Compute R2
  rss <- sum(r^2)
  resvar <- rss/rdf
  r.squared <- mss/(mss + rss)

# Copmute R2 with cell means codes ---------------------------------------------

  # Fit the lm with cell-means codes
  lm_cell <- lm(Petal.Length ~ Species - 1, data = iris)
  summary(lm_cell)

  # Extract basic object we need
  f   <- lm_cell$fitted.values # fitted values (y_hat)
  r   <- lm_cell$residuals     # residuals (y_true - y_hat)
  rdf <- lm_cell$df.residual   # residual degrees of freedom (n - k - 1)

  # By default, mean square error produced as if y had been centered
  mss <- sum(f^2)

  # We did not center y! So we want to subtract the mean of the fitted values
  mss <- sum((f - mean(f))^2)

  # And then we can compute the R2 as usual
  rss <- sum(r^2)
  resvar <- rss/rdf

  # Compute R2
  r.squared <- mss/(mss + rss)
