library(RcppArmadillo)
library(mvtnorm)

source("functions.R")
Rcpp::sourceCpp('nz.cpp')

sim <- simDiagModel(
  d = 5,
  n = 200,
  h2 = 0.5,
  r0 = 0.2,
  gamma0 = 0.0001,
  a. = 0.02,
  b. = 0.5,
  c. = 0.15,
  rho = 0.5,
  sigmaE2 = 0.01,
  nyear = 20,
  dt = 0.01,
  burnin = 10
)
