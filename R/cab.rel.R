#' The element $c_a$ or $c_b$ which is in the set S in relative case
#'
#' @param dist A characther string specifying the probability distribution. Distributions 
#' "binomial", "negative binomial" and "poisson".
#' @param ab A positive real number. The bound of the interval (a, b). Input 'a' or 'b'.
#' @param n A positive interger representing the sample size.
#' @param e A positive real number representing the maximum admissible relative estimation 
#' error.
#' @param w A positive real number representing the aliquot volume (only available for 
#' negative binomial distribution).
#' @param phi A positive real number representing the shape parameter for negative 
#' binomial distribution.
#'
#' @return A probability.
#'
#' @noRd
#'
cab.rel <- function(dist, ab, n, e, w, phi) {
  if (dist == "binomial") {
    out <- stats::pbinom(ceiling(n*ab*(1 + e)) - 1, size = n, 
                  prob = ab) - stats::pbinom(floor(n*ab*(1 - e)), size = n, prob = ab)
  }	
  if (dist == "negative binomial") {
    out <- stats::pnbinom(ceiling(n*w*ab*(1 + e)) - 1, size = n*phi, 
                   mu = n*w*ab) - stats::pnbinom(floor(n*w*ab*(1 - e)), size = n*phi, mu = n*w*ab)
  }
  if (dist == "poisson") {
		out <- stats::ppois(ceiling(n*ab*(1 + e)) - 1, n*ab) - stats::ppois(pmax(0, floor(n*ab*(1 - e))), n*ab)
	}	
	return(out)	
}
