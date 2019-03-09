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
  output <- matrix(NA, N, k)
  for (i in 1:N) {
    D <- stats::rgamma(k, shape = 1, rate = 1)
    output[i, ] <- D/sum(D)
  }
  return(output)
}
