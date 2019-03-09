#' Bayesian sample size in a decision-theoretic approach for the Poisson/Dirichlet process.
#'
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution.
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $\lam_0$ and shape parameter $\theta_0$.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param w A positive real number representing the aliquot volume.
#' @param c A positive real number representing the cost of colect one aliquot.
#' @param nmax A positive integer representing the maximum number for compute the Bayes risk.
#' Default is 100.
#' @param nrep A positive integer representing the number of samples taken for each $n$.
#' @param R Number of replicates used in the simulation. Default is 100.
#' @param plot Boolean. If TRUE (default) it plot the estimated Bayes risks and the fitted
#' curve.
#' @param ... Currently ignored. 
#'
#' @return An integer representing the sample size.
#' @export
#'
#' @examples
bss.dt.PoiDP <- function(lam0, theta0, alpha, w, c, nmax = 1E2, nrep = 1E1, R = 1E2, 
                     plot = TRUE, ...) {
  cl <- match.call()
  ns <- seq(3, nmax, by = 5)
  risk <- numeric()
    for (n in ns) {  
      for (i in 1:nrep) {
        loss <- numeric()
        for (j in 1:R) {
          x <- stats::rnbinom(n, mu = w*lam0, size = theta0)
          obj.vpost <- var.postDPmix(x = x, w = w, lam0 = lam0, theta0 = theta0, 
                                     alpha = alpha)
          loss <- append(loss, sum(obj.vpost$varp*dnorm(obj.vpost$lam, mean = 10, 
                                                        sd = 1E1)) + c*n)
        }
        risk <- append(risk, mean(loss))
      }
    }
  Y <- log(risk - c*rep(ns, each = nrep))
  mod <- stats::lm(Y ~ I(log(rep(ns + 1, each = nrep))))
  E <- as.numeric(exp(mod$coef[1]))
  G <- as.numeric(-mod$coef[2])
  nmin <- ceiling((E*G/c)^(1/(G + 1)) - 1)
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
