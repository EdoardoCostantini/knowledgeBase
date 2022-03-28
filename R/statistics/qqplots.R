### Title:    Stats & Methods Lab 6 QQ plots in depth
### Author:   Edoardo Costantini
### Created:  2021-03-25
### Tag:      MLR; Assumptions; Distributions

rm(list = ls())

library(MASS)     # For the 'Cars93' data

## Load some data:
data(Cars93)

## Fit some model
out1 <- lm(Price ~ Horsepower + MPG.city + Passengers, data = Cars93)
summary(out1)

## Save the residuals and fitted values from out1:
res1  <- resid(out1)
yHat1 <- predict(out1)

## Check the normality of residual via Q-Q Plot:
qqnorm(res1, ylim = c(-15, 35), xlim = c(-3, 3))
qqline(res1)

## Alternative using plot.lm function (resiuals are standardized)
plot(out1, which = 2)

## QQ plot logic

# Step 1 - Empirical Cumulative Propobability Distribution
  probs <- ppoints(res1)
  plot(x = sort(res1), y = probs)

# Step 2 - What quantiles do these values match on the normal PDF?
  theoQ <- qnorm(probs)

  # Theoretical Standard normal disitribution
  x <- seq(-4, 4, .01)
  plot(x, y = dnorm(x), type = "l")
    points(x = theoQ, y = dnorm(theoQ))
    abline(v = qnorm(probs))

# Step 3 - Create qqplot plot
  # Unstandardized residuals
  plot(x = theoQ, y = sort(res1),
       ylim = c(-15, 10), xlim = c(-3, 3))
  qqline(res1)

  # Standardized residuals
  plot(x = theoQ,
       y = sort(scale(res1)))
  qqline(scale(res1))

  plot(out1, which = 2)

  # What if residuals were actually normally distributed?
  # note the 45 degree angle
  set.seed(20220325)
  res_normal <- rnorm(1e4, 10, 2)
  theoQ <- qnorm(ppoints(res_normal))
  plot(x = theoQ,
       y = sort(res_normal))
  plot(x = theoQ,
       y = sort(scale(res_normal)))
