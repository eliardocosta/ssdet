#' Random generation for the Dirichlet distribution with all parameters equal to 1.
#'
#' @param N Number of samples to generate.
#' @param k Size of the samples.
#'
#' @return A matrix with the generated samples.
#' 
#' @noRd
#'
rdiri <- function(N, k) { 
  out.diri <- matrix(stats::rgamma(N*k, shape = 1, rate = 1 ), ncol = k, byrow = TRUE)
  return(out.diri/apply(out.diri, 1, sum))
}
