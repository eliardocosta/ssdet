#' Bayesian sample size in a decision-theoretic approach for the negative binomial/Pearson Type VI model
#'
#' @param lf 1 or 2, representing the loss function used.
#' @param lam0 A positive real number representing the prior expected value for the prior 
#' gamma distribution.
#' @param theta0 A positive real number representing the shape parameter for the prior 
#' gamma distribution.
#' @param phi A positive real number representing a scale parameter of the prior distribution.
#' @param w A positive real number representing the aliquot volume.
#' @param c A positive real number representing the cost of colect one aliquot.
#' @param rho A number in (0, 1). The probability of the credible interval is $1-rho$. Only
#' for lost function 1.
#' @param gam A positive real number connected with the credible interval when using lost
#' function 2.
#' @param nmax A positive integer representing the maximum number for compute the Bayes risk.
#' Default is 100.
#' @param nrep A positive integer representing the number of samples taken for each $n$.
#' @param lrep A positive integer representing the number of samples taken for $S_n$.
#' @param plot Boolean. If TRUE (default) it plot the estimated Bayes risks and the fitted
#' curve.
#' @param ... Currently ignored. 
#'
#' @return An integer representing the sample size.
#' @export
#'
#' @examples
bss.dt.nbPearson <- function(lf, lam0, theta0, phi, w, c, rho = NULL, gam = NULL, nmax = 1E2,
                             nrep = 1E1, lrep = 5E1, plot = TRUE, ...) {
  cl <- match.call()
  ns <- seq(1, nmax, by = 5)
  risk <- numeric()
  if (lf == 1) {
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        sn <- numeric()
        for (i in 1:lrep) {
          lam <- PearsonDS::rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, location = 0,
                                       scale = phi/w)
          x <- stats::rnbinom(length(lam), mu = w*lam, size = phi)
          sn <- append(sn, sum(x))
        }
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        a <- PearsonDS::qpearsonVI(rho/2, a = kappa, b = psi, location = 0, scale = phi/w)
        b <- PearsonDS::qpearsonVI(1 - rho/2, a = kappa, b = psi, location = 0, scale = phi/w)
        tau <- (b - a)/2
        loss <- (phi/w)*(
          PearsonDS::ppearsonVI(b, a = kappa + 1, b = psi - 1, location = 0, scale = phi/w,
                                lower.tail = FALSE) - PearsonDS::ppearsonVI(a, a = kappa + 1,
                                                                            b = psi - 1,
                                                                            location = 0,
                                                                            scale = phi/w)) + c*n
        risk <- append(risk, mean(loss))
      }
    }
  } else if (lf == 2){
    for (n in ns) {
      for (i in 1:nrep) {
        loss <- numeric()
        sn <- numeric()
        for (i in 1:lrep) {
          lam <- PearsonDS::rpearsonVI(n, a = theta0, b = theta0/lam0 + 1, 
                                       location = 0, scale = phi/w)
          x <- stats::rnbinom(n, mu = w*lam, size = phi)
          sn <- append(sn, sum(x))
        }
        kappa <- theta0 + sn
        psi <- theta0/lam0 + n*phi + 1
        qcon <- 1/(psi - 1)
        medpos <- (phi/w)*kappa/(psi - 1)
        varpos <- (kappa/(psi - 1)^2)*((medpos + 1)/(1 - qcon))*(phi/w)^2
        loss <- 2*sqrt(gam*varpos) + c*n
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
