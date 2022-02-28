# Directory description
Directory hosting an EM algorithm to estimate means and covariance matrices in the presence
of a general missing data pattern.

# Repository contents
This directory contains:
- **emSchafer.R** containing a function to run an EM algorithm coded as closely as possible to
  Schafer 1997 figure 5.2 (p.218)
- **emFast.R** containing a function to run an EM algorithm coded in a more R way.
- **em.R** a script running `emSchafer()` and `emFast()` on data examples from Little & Rubin 20--
  (multivariate monotone pattern) and Schafer 1997 (multivariate general missing data pattern).
- **sweepGoodnight.R** containing a function implementing a naive sweep operator
- **sweep.R** a script exploring the depths of the sweep operator in abstract 
  (no missing values involved)