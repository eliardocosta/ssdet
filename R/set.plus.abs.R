#' The element $c_{+}$ which is in the set S in absolute case
#'
#' @param dist A characther string specifying the probability distribution. Distributions 
#' "binomial", "negative binomial" and "poisson".
#' @param a A positive real number representing the lower bound of the interval (a, b).
#' @param b A positive real number representing the upper bound of the interval (a, b).
#' @param n A positive interger representing the sample size.
#' @param e A positive real number representing the maximum admissible absolute estimation 
#' error.
#' @param w A positive real number representing the aliquot volume (only available for 
#' negative binomial distribution).
#' @param phi A positive real number representing the shape parameter for negative 
#' binomial distribution.
#'
#' @return A set of probabilities.
#' 
#' @noRd
#'
set.plus.abs <- function(dist, a, b, n, e, w, phi) {
  if (dist == "binomial" || dist == "poisson") {
    lim1 <- max(0, 1 + floor(n*(a - e)))
    lim2 <- ceiling(n*(b - e)) - 1
  }
  if (dist == "negative binomial") {
		lim1 <- max(0, 1 + floor(n*w*(a - e)))
    lim2 <- ceiling(n*w*(b - e)) - 1
	}
	l <- seq.int(from = lim1, by = 1, length.out = lim2 - lim1 + 1)
	if (dist == "binomial") {
	  set <- stats::pbinom(l - 1 + ceiling(2*n*e), size = n, 
	                prob = l/n + e) - stats::pbinom(l, size = n, prob = l/n + e)
	}
	if (dist == "negative binomial") {
		set <- stats::pnbinom(l - 1 + ceiling(2*n*w*e), size = n*phi, 
		               mu = l + n*w*e) - stats::pnbinom(l, size = n*phi, mu = l + n*w*e)	
	}
	if (dist == "poisson") {
		set <- stats::ppois(l - 1 + ceiling(2*n*e), l + n*e) - stats::ppois(pmax(0, l), l + n*e)
	}
	return(set)
}
