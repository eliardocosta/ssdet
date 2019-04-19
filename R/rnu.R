#' Generate samples of the mixing distribution of the mixture of Dirichlet processes.
#'
#' @param nsam Number of samples to generate.
#' @param x Observed counts of the Poisson distribution of the model.
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
rnu <- function(nsam, x, w, lam0, theta0, alpha, nburn = 5E2) {
  n <- length(x)
  lam <- round(stats::rgamma(n, shape = theta0, rate = theta0/lam0), 3)
  out.nu <- matrix(nrow = nsam, ncol = n)
  for (t in seq_len(nburn)) {
    lam <- sapply(seq_len(n), function(i) {
      q0 <- alpha*stats::dnbinom(x[i], mu = w*lam0, size = theta0)
      qk <- stats::dpois(x[i], lambda = w*lam[-i])
      cn <- q0 + sum(qk) # normalization constant
      q0n <- q0/cn
      qkn <- qk/cn
      u <- stats::runif(1)
      if (u <= q0n) {
        lam[i] <- stats::rgamma(1, shape = theta0 + x[i], rate = (w + theta0/lam0))
      } else {
        lam[i] <- sample(x = lam[-i], size = 1, prob = qkn)
      }
      return(lam[i])
    })
    lam.table <- as.data.frame(table(lam), stringsAsFactors = FALSE)
    lam.star <- as.vector(lam.table[['lam']], mode = "numeric")
    n.star <- lam.table[['Freq']] 
    S.set <- numeric(n)
    S.set <- sapply(seq_len(n), function(i) which(lam[i] == lam.star))
    for (i in seq_len(length(n.star))) {
      lam[which(S.set == i)] <- stats::rgamma(1, shape = theta0 + sum(x[which(S.set == i)]), 
                                              rate = n.star[i] + theta0/lam0)
    }
  }
  for (t in seq_len(nsam)) {
    lam <- sapply(seq_len(n), function(i) {
      q0 <- alpha*stats::dnbinom(x[i], mu = w*lam0, size = theta0)
      qk <- stats::dpois(x[i], lambda = w*lam[-i])
      cn <- q0 + sum(qk) # normalization constant
      q0n <- q0/cn
      qkn <- qk/cn
      u <- stats::runif(1)
      if (u <= q0n) {
        lam[i] <- stats::rgamma(1, shape = theta0 + x[i], rate = (w + theta0/lam0))
      } else {
        lam[i] <- sample(x = lam[-i], size = 1, prob = qkn)
      }
      return(lam[i])
    })
    lam.table <- as.data.frame(table(lam), stringsAsFactors = FALSE)
    lam.star <- as.vector(lam.table[['lam']], mode = "numeric")
    n.star <- lam.table[['Freq']] 
    S.set <- numeric(n)
    S.set <- sapply(seq_len(n), function(i) which(lam[i] == lam.star))
    for (i in seq_len(length(n.star))) {
      lam[which(S.set == i)] <- stats::rgamma(1, shape = theta0 + sum(x[which(S.set == i)]), 
                                              rate = n.star[i] + theta0/lam0)
    }
    out.nu[t, ] <- lam
  }
  return(out.nu)
}
