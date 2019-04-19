#' Algorithm of Moller to simulate the functional of a Dirichlet process
#'
#' @param llam 
#' @param ulam 
#' @param eps 
#'
#' @return a simulated value.
#' @export
#'
#' @noRd
#' 
lam.moller.alg <- function(llam, ulam, eps, alpha, lam0, theta0) {
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
