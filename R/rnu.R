#' Generate samples of the mixing distribution of the mixture of Dirichlet processes.
#'
#' @param nsam Number of samples to generate.
#' @param x Counts of the Poisson distribution of the model.
#' @param w A positive real number representing the aliquot volume.
#' @param lam0 A positive real number representing a hyperparameter of the $F_0$ base distribution. 
#' @param theta0 A positive real number representing a hyperparameter of the $F_0$ base distribution. We consider $F_0$ as the gamma distribution with mean $\lam_0$ and shape parameter $\theta_0$.
#' @param alpha Shape parameter of the Dirichlet process.
#' @param nburn Number of burn-in samples.
#'
#' @return A matrix with samples of the mixing distribution of the mixture of Dirichlet processes.
#'
#' @noRd
#' 
rnu <- function(nsam, x, w, lam0, theta0, alpha, nburn = 1E3) {
  n <- length(x)
  lam <- round(rgamma(n, shape = theta0, rate = theta0/lam0), 3)
  output <- matrix(NA, nsam, n)
  nclu <- length(lam)
  nite <- nsam + nburn
  for (t in 1:nite) {
    for (i in 1:n) { 
      q0 <- alpha*dnbinom(x[i], mu = w*lam0, size = theta0)
      qk <- dpois(x[i], lambda = w*lam[-i])
      cn <- q0 + sum(qk) # constante de normalizacao
      q0n <- q0/cn
      qkn <- qk/cn
      u <- runif(1)
      lam[i] <- ifelse(u <= q0n, rgamma(1, shape = theta0 + x[i], rate = (w + theta0/lam0)),
                       sample(x = lam[-i], size = 1, prob = qkn))
    }
    lam <- round(lam, 3)
    if (t > nburn) {
      output[t - nburn, ] <- lam
    } 
    lam.table <- as.data.frame(table(lam), stringsAsFactors = FALSE)
    lam.star <- as.vector(lam.table[, 1], mode = "numeric")
    n.star <- as.vector(lam.table[, 2])
    S <- numeric()
    for (i in 1:n) {
      S <- append(S, which(lam[i] == lam.star))
    }
    for (i in 1:length(n.star)) {
      lam[which(S == i)] <- rgamma(1, shape = theta0 + sum(x[which(S == i)]), rate = n.star[i] + theta0/lam0)
    }
  }
  return(output)
}
