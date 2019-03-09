#' Expected posterior of the Poisson/mixture of Dirichlet process with a gamma distribution as the $F_0$ base distribution.
#'
#' @param x Observed counts of the Poisson distribution of the model.
#' @param w A positive real number representing the aliquot volume.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $\lam_0$ and shape parameter $\theta_0$.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param cgrid Length of the grid used in the points for which the probability is computed.
#' @param nsam Number of samples to use in the simulation.
#'
#' @return A list with cumulative posterior probabilities.
#'
#' @noRd
exp.postDPmix <- function(x = x, w = w, lam0 = lam0, theta0 = theta0, alpha = alpha, 
                          cgrid = 1E-2, nsam = 1E2) {
  n <- length(x)
  samcon.lam <- rnu(nsam = nsam, x = x, w = w, lam0 = lam0, theta0 = theta0, alpha = alpha)
  grid <- seq(0, ceiling(max(samcon.lam)), cgrid)
  probs <- matrix(NA, nsam, length(grid))
  probs[ ,1] <- rep(0, nsam)
  for (i in 1:nsam) {
    for (j in 2:length(grid)) {
      probs[i, j] <- (alpha/(alpha + n))*stats::pgamma(grid[j], shape = theta0, rate = theta0/lam0)+
        (1/(alpha + n))*(length(which(samcon.lam[i, ] <= grid[j])))
    }
  }
  return(list(lam = grid, cumprob = apply(probs, 2, mean), cgrid = cgrid))
} 
