#' Bayesian sample size using ACC or ALC criterion for the Poisson/gamma model
#'
#' @param crit A characther string specifying the criterion. Criteria: "ACC", "ALC" and
#' "ALCApprox"
#' @param lam0 A positive real number representing the prior expected value for the prior 
#' gamma distribution
#' @param theta0 A positive real number representing the shape parameter for the prior 
#' gamma distribution
#' @param w A positive real number representing the aliquot volume 
#' @param rho A number in (0, 1). The probability of the credible interval is equal or
#' greater than $1-rho$ depending on the criterion used.
#' @param len A positive real number representing the length of the credible intervals in the
#' ACC criterion.
#' @param len.max A positive real number representing the maximum length of the credible 
#' intervals in the ALC criterion. 
#' @param R Number of replicates used in the simulation. Default is 1000.
#' @param n0 A positive integer representing the initial sample size in which the function 
#' will check the criterion. Default is 1.
#'
#' @return An integer representing the sample size
#' @export
#'
bss.ci.poigam <- function(crit, lam0, theta0, w, rho, len = NULL, len.max = NULL, 
                          R = 1E3, n0 = 1) {
  cl <- match.call()
  if (crit == "ALCAprox") {
    zrho <- stats::qnorm(1 - rho/2)
    n <- (theta0/(w*lam0))*(((lam0/theta0)*(2*zrho/len.max)*(base::gamma(theta0 + 
                                                                     0.5)/base::gamma(theta0)))^2 - 1)
  }
  if (crit == "ACC") {
    cov <- 0 
    n <- n0
    while (mean(cov) < 1 - rho) {
      n <- n + 1
      cov <- numeric()
      for (i in 1:R) {
        sn <- stats::rnbinom(1, mu = n*w*lam0, size = n*theta0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        ab <- hpd.gamma(kappa = kappa, psi = psi, len = len)
        cov <- append(cov, stats::pgamma(ab[2], shape = kappa, 
                                  rate = psi) - stats::pgamma(ab[1], shape = kappa, rate = psi))
      }
    }
  } 
  if (crit == "ALC") { 
    len <- len.max + 1
    n <- n0
    while (base::mean(len) > len.max) {
      n <- n + 1
      len <- numeric()
      for (i in 1:R) {
        sn <- stats::rnbinom(1, mu = n*w*lam0, size = n*theta0)
        kappa <- theta0 + sn
        psi <- n*w + theta0/lam0
        ab <- hpd.gamma(kappa = kappa, psi = psi, rho = rho)
        len <- append(len, ab[2] - ab[1])
      }
    }
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nSample size:\n")
  cat("n  = ", n)
} 
