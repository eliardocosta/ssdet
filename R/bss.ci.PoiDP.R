#' Bayesian sample size using ACC or ALC criterion for the Poisson/Dirichlet process model
#'
#' @param crit A characther string specifying the criterion. Criteria: "ACC" and "ALC".
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution.
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $lam_0$ and shape parameter $theta_0$.
#' @param w A positive real number representing the aliquot volume.
#' @param rho A number in (0, 1). The probability of the credible interval is equal or
#' greater than $1-rho$ depending on the criterion used.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param len A positive real number representing the length of the credible intervals in the
#' ACC criterion.
#' @param len.max A positive real number representing the maximum length of the credible 
#' intervals in the ALC criterion. 
#' @param eps ...
#' @param cgrid Length of the grid used in the points for which the probability is computed.
#' @param R Number of replicates used in the simulation. Default is 100.
#' @param n0 A positive integer representing the initial sample size in which the function 
#' will check the criterion. Default is 1.
#' @param inc ...
#'
#' @return An integer representing the sample size.
#' @export
#'
bss.ci.PoiDP <- function(crit, lam0, theta0, w, rho, alpha, len = NULL, len.max = NULL, 
                     eps = NULL, cgrid = 1E-2, R = 1E2, n0 = 2, inc = c(1E2, 1E1, 5)) {
  cl <- match.call()
  if (crit == "ACC") { 
    n <- n0
    for (i in 1:length(inc)) {
      cov <- 0
      while (mean(cov) < 1 - rho) {
        n <- n + inc[i]
        cov <- numeric()
        for (j in 1:R) {
          x <- rnbinom(n, mu = w*lam0, size = theta0)
          obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, theta0 = theta0, 
                                     alpha = alpha, cgrid = cgrid)
          obj.dpost <- extract.prob(obj.ppost)
          vals <- obj.dpost$lam
          probs <- obj.dpost$ddist
          conj.vals <- numeric()
          conj.ind <- numeric()
          for (k in 1:ceiling(len/obj.ppost$cgrid)) {
            ind.max <- which.max(probs) 
            conj.vals <- append(conj.vals, vals[ind.max])
            conj.ind <- append(conj.ind, which(obj.dpost$lam == vals[ind.max]))
            vals <- vals[-ind.max]
            probs <- probs[-ind.max]
          }
          cov <- append(cov, sum(obj.dpost$ddist[conj.ind]))
        }
      }
      if (i < length(inc)) n <- n - inc[i]
    }
  } 
  if (crit == "ALC") {
    n <- n0
    for (i in 1:length(inc)) {
      len <- len.max + 1
      while (mean(len) > len.max) {
        n <- n + inc[i]
        len <- numeric()
        for (j in 1:R) {
          x <- stats::rnbinom(n, mu = w*lam0, size = theta0)
          obj.ppost <- exp.postDPmix(x = x, w = w, lam0 = lam0, theta0 = theta0, 
                                     alpha = alpha, cgrid = cgrid)
          obj.dpost <- extract.prob(obj.ppost)
          vals <- obj.dpost$lam
          probs <- obj.dpost$ddist
          conj.vals <- numeric()
          cov <- 0
          while (sum(cov) < 1 - rho) {
            ind.max <- which.max(probs) 
            conj.vals <- append(conj.vals, vals[ind.max])
            for (k in 1:length(conj.vals)) {
              cov[k] <- obj.dpost$ddist[obj.dpost$lam == conj.vals[k]]
            }
            vals <- vals[-ind.max]
            probs <- probs[-ind.max]
          }
          len <- append(len, length(conj.vals)*obj.ppost$cgrid)
        }
      }
      if (i < length(inc)) n <- n - inc[i]
    }
  } 
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nSample size:\n")
  cat("n  = ", n, "\n")
} 
