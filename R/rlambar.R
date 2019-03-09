#' Random generation for the functional mean of the Dirichlet process with a gamma distribution as the $F_0$ base distribution.
#'
#' @param N Number of observations.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $lam_0$ and shape parameter $theta_0$.
#' @param eps Tolerance limit used in the simulation algorithm.
#'
#' @return A random sample of the functional mean of the Dirichlet process.
#' @export
#'
rlambar <- function(N, alpha, lam0, theta0, eps = 5E-2) {
  output <- numeric()
  for (i in 1:N) {
    B <- stats::rbeta(1, 1, alpha)
    xi <- stats::rgamma(1, shape = theta0, rate = theta0/lam0)
    ulam <- B*xi + (1 - B)*.Machine$double.xmax
    llam <- B*xi
    j <- 1
    while (abs(ulam - llam) > eps) {
      B <- stats::rbeta(1, 1, alpha)
      xi <- stats::rgamma(1, shape = theta0, rate = theta0/lam0)
      ulam <- B*xi + (1 - B)*ulam
      llam <- B*xi + (1 - B)*llam
      j <- j + 1
    }
    output <- append(output, ulam)
  }
  return(output)
} 
