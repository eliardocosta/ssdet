#' Random generation for the Dirichlet process with a gamma distribution as the $F_0$ base distribution.
#'
#' @param n Sample size.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution.
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $lam_0$ and shape parameter $theta_0$.
#' @param eps Tolerance limit used in the simulation algorithm. Default is 0.01.
#'
#' @return A sample of size n from the Dirichlet process.
#' @export
#'
rlam <- function(n, alpha, lam0, theta0, eps = 1E-2) {
  M <- ceiling(1 - alpha*log(eps/(4*n))) # Ishwaran & James (2001, Theo. 2)
  lam <- numeric()
  p <- numeric()
  V <- c(stats::rbeta(M - 1, 1, alpha), 1)
  p[1] <- V[1]
  for (i in 2:M) {
    p[i] <- V[i]*(1 - V[i - 1])*p[i - 1]/V[i - 1]
  }
  sam.lam <- stats::rgamma(M, shape = theta0, rate = theta0/lam0)
  lam <- sample(sam.lam, n, replace = TRUE, prob = p)
  return(lam)
}
