#' Compute the HPD interval for the posterior distribution of the negative binomial/Pearson Type VI model.
#'
#' @param kappa Fist parameter of the posterior distribution.
#' @param psi Second parameter of the posterior distribution.
#' @param phi A positive real number representing a scale parameter of the prior distribution.
#' @param w A positive real number representing the aliquot volume.
#' @param len Length of the HPD interval
#' @param rho Probability of the HPD interval
#'
#' @return
#'
#' @noRd
#' 
hpdPearsonVI <- function(kappa, psi, phi, w, len = NULL, rho = NULL) {
  if (is.null(len)) {
    fun <- function(x) c(F1 = PearsonDS::ppearsonVI(x[2], a = kappa, b = psi, location = 0,
                        scale = phi/w) - PearsonDS::ppearsonVI(x[1], a = kappa, b = psi,
                        location = 0, scale = phi/w) - 1 + rho,
                        F2 = PearsonDS::dpearsonVI(x[1], a = kappa, b = psi, location = 0,
                        scale = phi/w) - PearsonDS::dpearsonVI(x[2], a = kappa, b = psi,
                        location = 0, scale = phi/w))
    roots <- c(-2, -1)
    starts <- PearsonDS::qpearsonVI(c(rho/2, 1 - rho/2), a = kappa, b = psi, location = 0,
                                    scale = phi/w)
    while (all(roots <= 0)) {
      sol.finder <- rootSolve::multiroot(f = fun, start = starts)
      roots <- sol.finder$root
      starts <- starts + 0.1
    }
  }
  if (is.null(rho)) {
    fun <- function(x) {
      (kappa - 1)*log(1 + len/x) - (kappa + psi)*log(1 + w*len/(phi + w*x))
    }
    sol <- uniroot(fun, lower = 1E-1, upper = 1E2, extendInt = "downX")
    roots <- c(sol$root, sol$root + len)
  }
  return(roots)
} # FIM
