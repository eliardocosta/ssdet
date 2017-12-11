#' The element $c_{-}$ which is in the set S in relative case
#'
#' @param dist A characther string specifying the probability distribution. Distributions 
#' "binomial", "negative binomial" and "poisson"
#' @param a A positive real number representing the lower bound of the interval (a, b)
#' @param b A positive real number representing the upper bound of the interval (a, b)
#' @param n A positive interger representing the sample size
#' @param e A positive real number representing the maximum admissible relative estimation 
#' error
#' @param w A positive real number representing the aliquot volume (only available for 
#' negative binomial distribution)
#' @param phi A positive real number representing the shape parameter for negative 
#' binomial distribution
#'
#' @return A set of probabilities
#' 
#' @noRd
#'
set.minus.rel <- function(dist, a, b, n, e, w, phi) {
  if (dist == "binomial" || dist == "poisson") {
    lim1 <- 1 + floor(n*a*(1 + e))
    lim2 <- ceiling(n*b*(1 + e)) - 1
  }
  if (dist == "negative binomial") {
		lim1 <- 1 + floor(n*w*a*(1 + e))
    lim2 <- ceiling(n*w*b*(1 + e)) - 1
	}
	l <- seq.int(from = lim1, by = 1, length.out = lim2 - lim1 + 1)
	if (dist == "binomial") {
	  set <- stats::pbinom(l - 1, size = n, prob = l/(n*(1 + e))) - stats::pbinom(floor(l*(1 - e)/(1 + e)), 
	                                                                size = n, prob = l/(n*(1 + e)))
	}
	if (dist == "negative binomial") {
		set <- stats::pnbinom(l - 1, size = n*phi, 
		               mu = l/(1 + e)) - stats::pnbinom(floor(l*(1 - e)/(1 + e)), 
		                                         size = n*phi, mu = l/(1 + e))
	}
	if (dist == "poisson") {
		set <- stats::ppois(l - 1, l/(1 + e)) - stats::ppois(pmax(0, floor(l*(1 - e)/(1 + e))), l/(1 + e))
	}
	return(set)
}
