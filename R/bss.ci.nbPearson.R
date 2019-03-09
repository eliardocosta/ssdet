#' Bayesian sample size using ACC or ALC criterion for the negative binomial/Pearson Type VI model
#'
#' @param crit A characther string specifying the criterion. Criteria: "ACC" and "ALC".
#' @param lam0 A positive real number representing a shape parameter of the prior distribution.
#' @param theta0 A positive real number representing a shape parameter of the prior distribution.
#' @param phi A positive real number representing a scale parameter of the prior distribution.
#' @param w A positive real number representing the aliquot volume. 
#' @param rho A number in (0, 1). The probability of the credible interval is equal or greater than $1-rho$ depending on the criterion used.
#' @param len A positive real number representing the length of the credible intervals in the ACC criterion.
#' @param len.max A positive real number representing the maximum length of the credible  intervals in the ALC criterion.
#' @param R Number of replicates used in the simulation. Default is 1000.
#' @param n0 A positive integer representing the initial sample size in which the function will check the criterion. Default is 1.
#'
#' @return An integer representing the sample size.
#' @export
#'
#' @examples
bss.ci.nbPearson <- function(crit, lam0, theta0, phi, w, rho, len = NULL, len.max = NULL, 
                             R = 1E3, n0 = 1) {
  cl <- match.call()
  if (crit == "ACC") {
    cov <- 0 
    n <- n0
    while (mean(cov) < 1 - rho) {
      n <- n + 1
      cov <- numeric()
      probs <- numeric()
      for (i in 1:R) {
        lam <- PearsonDS::rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, 
                                     scale = phi/w)
        x <- stats::rnbinom(n, mu = w*lam, size = phi)
        sn <- sum(x)
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        ab <- hpd.PearsonVI(kappa = kappa, psi = psi, phi = phi, w = w, len = len)
        cov <- append(cov, 
                      PearsonDS::ppearsonVI(ab[2], a = kappa, b = psi, location = 0, 
                                            scale = phi/w) - PearsonDS::ppearsonVI(ab[1], 
                                            a = kappa, b = psi, location = 0, scale = phi/w))
      }
    }
  } 
  if (crit == "ALC") {
    len <- len.max + 1
    n <- n0
    while (mean(len) > len.max) {
      n <- n + 1
      len <- numeric()
      for (i in 1:R) {
        lam <- PearsonDS::rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0, 
                                     scale = phi/w)
        x <- stats::rnbinom(n, mu = w*lam, size = phi)
        sn <- sum(x)
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        ab <- hpd.PearsonVI(kappa = kappa, psi = psi, phi = phi, w = w, rho = rho)
        len <- append(len, ab[2] - ab[1])
      }
    }
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nSample size:\n")
  cat("n  = ", n, "\n")
} 
