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
stats::qqnorm
qqline(res1)

## QQ plot logic

# Step 1 - Give each point its own quantile
probs <- ppoints(res1)

# Step 2 - What quantiles do this match on the probability distirbution?
theoQ <- qnorm(probs)
  # Theoretical Standard normal disitribution
  x <- seq(-4, 4, .01)
  plot(x, y = dnorm(x), type = "l")
  abline(v = qnorm(probs))

# Step 3 - Cross plot
plot(x = theoQ, y = sort(res1),
     ylim = c(-15, 35), xlim = c(-3, 3))
qqnorm(res1,
       ylim = c(-15, 35), xlim = c(-3, 3))
