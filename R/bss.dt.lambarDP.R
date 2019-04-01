#' Bayesian sample size in a decision-theoretic approach for the functional mean of the Dirichlet process with a gamma distribution as the $F_0$ base distribution.
#'
#' @param lf 1 or 2, representing the loss function used.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $lam_0$ and shape parameter $theta_0$.
#' @param w A positive real number representing the aliquot volume.
#' @param c A positive real number representing the cost of colect one aliquot.
#' @param rho A number in (0, 1). The probability of the credible interval is $1-rho$. Only
#' for lost function 1.
#' @param gam A positive real number connected with the credible interval when using lost
#' function 2.
#' @param nmax A positive integer representing the maximum number for compute the Bayes risk.
#' Default is 100.
#' @param nrep A positive integer representing the number of samples taken for each $n$.
#' @param lrep A positive integer representing the number of samples taken for $S_n$. Default is 50.
#' @param plot Boolean. If TRUE (default) it plot the estimated Bayes risks and the fitted
#' curve.
#' @param ... Currently ignored.
#'
#' @return An integer representing the sample size.
#' @export
#'
bss.dt.lambarDP <- function(lf, alpha, lam0, theta0, w, c, rho = NULL, gam = NULL, 
                            nmax = 1E2, nrep = 1E1, lrep = 5E1, plot = FALSE, ...) {
  cl <- match.call()
  ns <- seq(3, nmax, by = 5)
  risk <- numeric()
  if (lf == 1) {
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        for (j in 1:lrep) {
          x <- stats::rnbinom(n, mu = w*lam0, size = theta0)
          lam.xn <- rlambar.xn(N = 1E2, alpha = alpha, x = x, w = w, lam0 = lam0, 
                               theta0 = theta0)
          qs <- stats::quantile(lam.xn, probs = c(rho/2, 1 - rho/2))
          loss <- append(loss, sum(lam.xn[which(lam.xn > qs[2])])/1E2 - sum(lam.xn[which(lam.xn < qs[1])])/1E2 + c*n) 
        }
        risk <- append(risk, mean(loss))
      }
    }
  } else if (lf == 2) {
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        for (j in 1:lrep) {
          x <- stats::rnbinom(n, mu = w*lam0, size = theta0)
          lam.xn <- rlambar.xn(N = 5E1, alpha = alpha, x = x, w = w, lam0 = lam0, 
                               theta0 = theta0)
          loss <- append(loss, 2*sqrt(gam*stats::var(lam.xn)) + c*n)
        }
        risk <- append(risk, mean(loss))
      }
    }
  }
  Y <- log(risk - c*rep(ns, each = nrep))
  mod <- stats::lm(Y ~ I(log(rep(ns + 1, each = nrep))))
  E <- as.numeric(exp(mod$coef[1]))
  G <- as.numeric(-mod$coef[2])
  nmin <- ceiling((E*G/c)^(1/(G + 1))-1)
  if (plot == TRUE) {
    plot(rep(ns, each = nrep), risk, xlim = c(0, nmax), xlab = "n", ylab = "TC(n)")
    curve <- function(x) {c*x + E/(1 + x)^G}
    plot(function(x)curve(x), 0, nmax, col = "blue", add = TRUE)
    graphics::abline(v = nmin, col = "red")
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nSample size:\n")
  cat("n  = ", nmin, "\n")
}
