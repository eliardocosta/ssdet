#' Random generation for the posterior functional mean of the Poisson/Dirichlet process with a gamma distribution as the $F_0$ base distribution.
#'
#' @param N Number of observations.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param x Observed counts of the Poisson distribution of the model.
#' @param w A positive real number representing the aliquot volume.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $\lam_0$ and shape parameter $\theta_0$.
#'
#' @return A random sample of the posterior functional mean of the Poisson/Dirichlet process.
#' @export
#'
#' @examples 
rlambar.xn <- function(N, alpha, x, w, lam0, theta0) {
  Z <- rnu(nsam = N, x = x, w = w, lam0 = lam0, theta0 = theta0, alpha = alpha)
  D <- rdiri(N = N, k = length(x))
  lam <- rlambar(N = N, alpha = alpha, lam0 = lam0, theta0 = theta0)
  B <- rbeta(N, alpha, length(x))
  output <- B*lam + (1 - B)*apply(D*Z, 1, sum)
  return(output)
}
