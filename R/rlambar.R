#' Random generation for the functional mean of the Dirichlet process with a gamma distribution as the $F_0$ base distribution.
#'
#' @param N Number of observations.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $lam_0$ and shape parameter $theta_0$.
#' @param eps Tolerance limit used in the simulation algorithm.
#' @param ncore Number of cores to use in parallel computin. If NULL the function uses 1 core if there is only one core, if there is more than one cores uses one half of the cores.
#'
#' @return A random sample of the functional mean of the Dirichlet process.
#' @importFrom foreach "%dopar%"
#' @export
#'
rlambar <- function(N, alpha, lam0, theta0, ncore, eps = 1E-1) {
  B <- stats::rbeta(1, 1, alpha)
  xi <- stats::rgamma(1, shape = theta0, rate = theta0/lam0)
  ulam <- B*xi + (1 - B)*.Machine$double.xmax
  llam <- B*xi
  if (is.null(ncore)) {
    num.cores <- ceiling(parallel::detectCores()/2)
  } else {
    num.cores <- ncore
  }
  doParallel::registerDoParallel(num.cores)
  out.lambar <- foreach::foreach (i = 1:N, .combine = 'c') %dopar% {
    lambar.moller.alg(llam = llam, ulam = ulam, eps = eps, alpha = alpha, lam0 = lam0, 
                   theta0 = theta0)
    }
  doParallel::stopImplicitCluster()
  return(out.lambar)
}
