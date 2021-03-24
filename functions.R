createDiagModelMatrices <- function(a., b., c., h2, d) {
  # Appendix F
  I <- diag(d)
  P <- I
  G <- h2 * P
  p <- solve(P)
  
  a <- a. * p
  b <- b. * p
  c <- c. * p
  
  q. <- (1 + c. + b.)
  v. <- (1 + c. - c. ^ 2 / q.)
  J. <- b. ^ 2 / q. + b. ^ 2 * c. ^ 2 / (q. ^ 2 * v.)
  K. <- c. ^ 2 / q. + (c. - c. ^ 2 / q.) ^ 2 / v.
  L. <- b. * c. / q. - b. * c. * (c. - c. ^ 2 / q.) / (q. * v.)
  
  q <- q. * I
  v <- v. * I
  J <- J. * I
  K <- K. * I
  L <- L. * I
  DET <- (q. * v.) ^ (-d / 2)
  
  A <- solve(a)
  
  bJL <- b - J - L
  cKL <- c - K - L
  
  half_sum_elem_prod_P_a <- 0.5 * sum(P * a)
  
  ret <-
    list(
      P = P,
      G = G,
      a = a,
      b = b,
      c = c,
      J = J,
      K = K,
      L = L,
      DET = DET,
      A = A,
      bJL = bJL,
      cKL = cKL,
      half_sum_elem_prod_P_a = half_sum_elem_prod_P_a
    )
  ret
}

simDiagModel <-
  function(d,
           n,
           h2,
           r0,
           gamma0,
           a.,
           b.,
           c.,
           rho,
           sigmaE2,
           nyear,
           dt,
           burnin) {
    # nyear : total number of years
    # dt : continous approximation - divide each year into 1/dt time steps
    # burnin : number of years without extinction barrier
    if (!(1 / dt == round(1 / dt))) {
      stop("The number of time steps per year (1/dt) should be an integer, set e.g. dt = 0.01")
    }
    
    nu <- mu <- rep(0, d)
    cMats <- createDiagModelMatrices(a., b., c., h2, d)
    c0Mats <- createDiagModelMatrices(a., b., c. = 0, h2, d)
    
    N <- rep(r0 / gamma0 / n, n)
    z <- t(rmvnorm(n, mu, cMats$A * 0.1))
    dt <- 0.01
    sigmaE <- sqrt(sigmaE2)
    half_sigmaE2 <- 0.5 * sigmaE2
    
    zArr <- array(NA, dim = c(nyear, n, d))
    Nmat <- matrix(NA, nyear, n)
    ns <- matrix(NA, nyear, 2)
    teller <- 0
    indx <- 1:n
    time <- rep(NA, nyear)
    tid <- 0
    attach(cMats, pos = -1)
   
    
    for (i in 1:burnin) {
      res <-
        calcGamma(
          nu,
          z,
          b,
          c,
          J,
          c0Mats$J,
          K,
          L,
          gamma0,
          DET,
          c0Mats$DET,
          a,
          mu,
          d,
          N,
          bJL,
          c0Mats$bJL,
          cKL,
          G,
          dt,
          half_sum_elem_prod_P_a,
          r0,
          half_sigmaE2,
          sigmaE,
          rho,
          niter = 1 / dt
        )
      N <- res$N
      z <- res$z
      Nmat[i, ] <- N
      zArr[i, , ] <- z
      cat(i, "\n")
      
    }
    
    for (i in (burnin + 1):nyear) {
      res <-
        calcGamma(
          nu,
          z,
          b,
          c,
          J,
          c0Mats$J,
          K,
          L,
          gamma0,
          DET,
          c0Mats$DET,
          a,
          mu,
          d,
          N,
          bJL,
          c0Mats$bJL,
          cKL,
          G,
          dt,
          half_sum_elem_prod_P_a,
          r0,
          half_sigmaE2,
          sigmaE,
          rho,
          niter = 1 / dt
        )
      N <- res$N
      z <- res$z
      notExt <- which(N > 0.5)
      indx <- indx[notExt]
      N <- N[notExt]
      z <- z[, notExt]
      zArr[i, indx, ] <- t(z)
      Nmat[i, indx] <- N
      cat(i, "\n")
      
    }
  
    return(list(Nmat = Nmat, zArr = zArr))
    
  }
