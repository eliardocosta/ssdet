#' Algorithm of Moller to simulate the functional of a Dirichlet process
#'
#' @param llam lower value in the Markov chain.
#' @param ulam upper value in the Markov chain.
#' @param eps Tolerance limit used in the simulation algorithm.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $lam_0$ and shape parameter $theta_0$.
#'
#' @return a simulated value of the functional of a Dirichlet process.
#'
#' @noRd
#' 
lambar.moller.alg <- function(llam, ulam, eps, alpha, lam0, theta0) {
  j <- 1
  while (abs(ulam - llam) > eps) {
    B <- stats::rbeta(1, 1, alpha)
    xi <- stats::rgamma(1, shape = theta0, rate = theta0/lam0)
    ulam <- B*xi + (1 - B)*ulam
    llam <- B*xi + (1 - B)*llam
    j <- j + 1
  }
  return(ulam)
}
