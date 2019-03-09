#' Compute the HPD interval for the posterior distribution of the Poisson/gamma model
#'
#' @param kappa Fist parameter of the posterior distribution.
#' @param psi Second parameter of the posterior distribution.
#' @param rho Probability of the HPD interval.
#' @param len Length of the HPD interval.
#'
#' @return Two number representing the HPD interval.
#' 
#' @noRd
#'
hpd.gamma <- function(kappa, psi, rho = NULL, len = NULL) {
  if (is.null(len)) {
    fun <- function(x) c(F1 = stats::pgamma(x[2], shape = kappa, rate = psi) - stats::pgamma(x[1], shape = kappa, rate = psi) - 1 + rho,
                         F2 = stats::dgamma(x[1], shape = kappa, rate = psi) - stats::dgamma(x[2], shape = kappa, rate = psi))
    roots <- c(-1, -1)
    starts <- stats::qgamma(c(rho/2, 1 - rho/2), shape = kappa, rate = psi)
    while (all(roots <= 0)) {
      sol.finder <- rootSolve::multiroot(f = fun, start = starts)
      roots <- sol.finder$root
      starts <- starts + 0.1
    }
  }
  if (is.null(rho)) {
    a <- len/(exp((psi*len)/(kappa - 1)) - 1)
    roots <- c(a, a +len)
  }
  return(roots)  
} 
