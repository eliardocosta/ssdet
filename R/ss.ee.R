#' Sample size for controling the estimation error controling the maximum absolute error,
#' the maximum relative error or both
#'
#' @param dist A characther string specifying the probability distribution. Distributions: 
#' "binomial", "negative binomial" and "poisson".
#' @param rho A number in (0, 1). The probability 1-rho represents the minimum 
#' confidence level.
#' @param a A positive real number representing the lower bound of the interval (a, b).
#' @param b A positive real number representing the upper bound of the interval (a, b).
#' @param ea A positive real number representing the maximum admissible absolute estimation 
#' error.
#' @param er A positive real number representing the maximum admissible relative estimation 
#' error.
#' @param w A positive real number representing the aliquot volume (only available for 
#' negative binomial distribution). Default is NULL.
#' @param phi A positive real number representing the shape parameter for negative 
#' binomial distribution. Default is NULL.
#' @param n0 A positive integer representing the initial sample size in which the function 
#' will check the criterion. Default is 2.
#' @param lag A positive integer representing the increment which 'n' receives until reach
#' the criterion. Default is 1.
#'
#' @return An integer representing the sample size.
#' 
#' @references Costa, E. G., Lopes, R. M. & Singer, J. M. (2016). Sample size for estimating 
#' the mean concentration of organisms in ballast water. Journal of Environmental Management 
#' 180, 433–438.
#' @references Chen, X. (2011). Exact computation of minimum sample size for estimation of 
#' binomial parameters. Journal of Statistical Planning and Inference 141, 2622–2632.
#' @references Chen, Z. & Chen, X. (2016). Exact calculation of minimum sample size for 
#' estimating a Poisson parameter. Communications in Statistics-Theory and Methods 45, 
#' 4692–4715
#' 
#' @export
#'
ss.ee <- function(dist, rho, a, b, ea = NULL, er = NULL, w = NULL, phi = NULL, 
                  n0 = 2, lag = 1) {
  cl <- match.call()

  # Initial checks
  if (is.null(ea) && is.null(er)) {
    stop("At least one error must be specified")
  }
  if (any(rho < 0) || any(rho > 1)) {
    stop("The elements of vector 'rho' must lie in interval (0, 1)")
  }
  if (a >= b) {
    stop("'a' must be smaller than 'b'")
  }
  if (dist == "binomial") {
    if (any(a, b) < 0 || any(a, b) > 1) {
      stop("'a' and 'b' must lie in interval (0, 1)")
    }
  }
  if (dist == "negative binomial") {
    if (is.null(w) || is.null(phi)) {
      stop("'w' and 'phi' must be specified")
    }
  }
  rho <- sort(rho, decreasing = TRUE)
  n <- n0
  
  #  Computing mixed case (absolute and relative errors)
  if (!is.null(ea) && !is.null(er)) {
    if(any(c(ea, er) < 0) || er > 1) {
    stop("'ea' must be positive and 'er' must lie in interval (0, 1)")
    }
    if (ea/er <= a || ea/er >= b) {
      stop("The ratio 'ea/er' must lie in interval (a, b)")
    }
    cat("Calculating...\n")
    sizes <- numeric()
    for (i in 1:length(rho)) {
      setprobs <- c(cab.abs(dist = dist, ab = a, n = n, e = ea, w = w, phi = phi), 
                    cab.abs(dist = dist, ab = ea/er, n = n, e = ea, w = w, phi = phi), 
                    set.plus.abs(dist = dist, a = a, b = ea/er, n = n, e = ea, w = w, phi = phi),
                    set.minus.abs(dist = dist, a = a, b = ea/er, n = n, e = ea, w = w, phi = phi), 
                    cab.rel(dist = dist, ab = ea/er, n = n, e = er, w = w, phi = phi), 
                    cab.rel(dist = dist, ab = b, n = n, e = er, w = w, phi = phi), 
                    set.plus.rel(dist = dist, a = ea/er, b = b, n = n, e = er, w = w, phi = phi), 
                    set.minus.rel(dist = dist, a = ea/er, b = b, n = n, e = er, w = w, phi = phi)
                    )
      while (min(setprobs) <= 1 - rho[i]) {
        n <- n + lag
        setprobs <- c(cab.abs(dist = dist, ab = a, n = n, e = ea, w = w, phi = phi), 
                      cab.abs(dist = dist, ab = ea/er, n = n, e = ea, w = w, phi = phi), 
                      set.plus.abs(dist = dist, a = a, b = ea/er, n = n, e = ea, w = w, phi = phi),
                      set.minus.abs(dist = dist, a = a, b = ea/er, n = n, e = ea, w = w, phi = phi), 
                      cab.rel(dist = dist, ab = ea/er, n = n, e = er, w = w, phi = phi), 
                      cab.rel(dist = dist, ab = b, n = n, e = er, w = w, phi = phi), 
                      set.plus.rel(dist = dist, a = ea/er, b = b, n = n, e = er, w = w, phi = phi), 
                      set.minus.rel(dist = dist, a = ea/er, b = b, n = n, e = er, w = w, phi = phi)
                      )
      }
      sizes <- append(sizes, n)
    }
    if (lag > 1) {
      sizes <- sizes - lag
      for (i in 1:length(rho)) {
        setprobs <- c(cab.abs(dist = dist, ab = a, n = n, e = ea, w = w, phi = phi), 
                      cab.abs(dist = dist, ab = ea/er, n = n, e = ea, w = w, phi = phi), 
                      set.plus.abs(dist = dist, a = a, b = ea/er, n = n, e = ea, w = w, phi = phi),
                      set.minus.abs(dist = dist, a = a, b = ea/er, n = n, e = ea, w = w, phi = phi), 
                      cab.rel(dist = dist, ab = ea/er, n = n, e = er, w = w, phi = phi), 
                      cab.rel(dist = dist, ab = b, n = n, e = er, w = w, phi = phi), 
                      set.plus.rel(dist = dist, a = ea/er, b = b, n = n, e = er, w = w, phi = phi), 
                      set.minus.rel(dist = dist, a = ea/er, b = b, n = n, e = er, w = w, phi = phi)
                      )
        while (min(setprobs) <= 1 - rho[i]) {
          sizes[i] <- sizes[i] + 1
          setprobs <- c(cab.abs(dist = dist, ab = a, n = sizes[i], e = ea, w = w, phi = phi), 
                        cab.abs(dist = dist, ab = ea/er, n = sizes[i], e = ea, w = w, phi = phi), 
                        set.plus.abs(dist = dist, a = a, b = ea/er, n = sizes[i], e = ea, w = w, phi = phi),
                        set.minus.abs(dist = dist, a = a, b = ea/er, n = sizes[i], e = ea, w = w, phi = phi), 
                        cab.rel(dist = dist, ab = ea/er, n = sizes[i], e = er, w = w, phi = phi), 
                        cab.rel(dist = dist, ab = b, n = sizes[i], e = er, w = w, phi = phi), 
                        set.plus.rel(dist = dist, a = ea/er, b = b, n = sizes[i], e = er, w = w, phi = phi), 
                        set.minus.rel(dist = dist, a = ea/er, b = b, n = sizes[i], e = er, w = w, phi = phi)
                        )
        }
      } 
    }
  }
  
  # Computing absolute error case
  if (!is.null(ea) && is.null(er)) {
    if (ea < 0) {
      stop("The absolute error 'ea' must be positive")
    } 
    if (ea > b - a) {
      stop("'ea' must be smaller than the amplitude of the interval (a, b)")
    }
    cat("Calculating...\n")
    sizes <- numeric()
    for (i in 1:length(rho)) {
      setprobs <- c(cab.abs(dist = dist, ab = a, n = n, e = ea, w = w, phi = phi), 
                    cab.abs(dist = dist, ab = b, n = n, e = ea, w = w, phi = phi), 
                    set.plus.abs(dist = dist, a = a, b = b, n = n, e = ea, w = w, phi = phi), 
                    set.minus.abs(dist = dist, a = a, b = b, n = n, e = ea, w = w, phi = phi)
                    )
      while (min(setprobs) <= 1 - rho[i]) {
        n <- n + lag
        setprobs <- c(cab.abs(dist = dist, ab = a, n = n, e = ea, w = w, phi = phi), 
                      cab.abs(dist = dist, ab = b, n = n, e = ea, w = w, phi = phi), 
                      set.plus.abs(dist = dist, a = a, b = b, n = n, e = ea, w = w, phi = phi), 
                      set.minus.abs(dist = dist, a = a, b = b, n = n, e = ea, w = w, phi = phi)
                      )
      }
      sizes <- append(sizes, n)
    }
    if (lag > 1) {
      sizes <- sizes - lag
      for (i in 1:length(rho)) {
        setprobs <- c(cab.abs(dist = dist, ab = a, n = sizes[i], e = ea, w = w, phi = phi), 
                      cab.abs(dist = dist, ab = b, n = sizes[i], e = ea, w = w, phi = phi), 
                      set.plus.abs(dist = dist, a = a, b = b, n = sizes[i], e = ea, w = w, phi = phi), 
                      set.minus.abs(dist = dist, a = a, b = b, n = sizes[i], e = ea, w = w, phi = phi)
                      )
        while (min(setprobs) <= 1 - rho[i]) {
          sizes[i] <- sizes[i] + 1
          setprobs <- c(cab.abs(dist = dist, ab = a, n = sizes[i], e = ea, w = w, phi = phi), 
                        cab.abs(dist = dist, ab = b, n = sizes[i], e = ea, w = w, phi = phi), 
                        set.plus.abs(dist = dist, a = a, b = b, n = sizes[i], e = ea, w = w, phi = phi), 
                        set.minus.abs(dist = dist, a = a, b = b, n = sizes[i], e = ea, w = w, phi = phi)
                        )
        }
      }
    }
  }
  
  #  Computing relative error case
  if (is.null(ea) && !is.null(er)) {
    if (er < 0 || er > 1) {
      stop("The relative error 'er' must lie in interval (0, 1)")
    }
    cat("Calculating...\n")
    sizes <- numeric()
    for (i in 1:length(rho)) {
      setprobs <- c(cab.rel(dist = dist, ab = a, n = n, e = er, w = w, phi = phi), 
                    cab.rel(dist = dist, ab = b, n = n, e = er, w = w, phi = phi), 
                    set.plus.rel(dist = dist, a = a, b = b, n = n, e = er, w = w, phi = phi), 
                    set.minus.rel(dist = dist, a = a, b = b, n = n, e = er, w = w, phi = phi)
                    )
      while (min(setprobs) <= 1 - rho[i]) {
        n <- n + lag
        setprobs <- c(cab.rel(dist = dist, ab = a, n = n, e = er, w = w, phi = phi), 
                      cab.rel(dist = dist, ab = b, n = n, e = er, w = w, phi = phi), 
                      set.plus.rel(dist = dist, a = a, b = b, n = n, e = er, w = w, phi = phi), 
                      set.minus.rel(dist = dist, a = a, b = b, n = n, e = er, w = w, phi = phi)
                      )
      }
      sizes <- append(sizes, n)
    }
    if (lag > 1) {
      sizes <- sizes - lag
      for (i in 1:length(rho)) {
        setprobs <- c(cab.rel(dist = dist, ab = a, n = sizes[i], e = er, w = w, phi = phi), 
                      cab.rel(dist = dist, ab = b, n = sizes[i], e = er, w = w, phi = phi), 
                      set.plus.rel(dist = dist, a = a, b = b, n = sizes[i], e = er, w = w, phi = phi), 
                      set.minus.rel(dist = dist, a = a, b = b, n = sizes[i], e = er, w = w, phi = phi)
                      )
        while(min(setprobs) <= 1 - rho[i]) {
          sizes[i] <- sizes[i] + 1
          setprobs <- c(cab.rel(dist = dist, ab = a, n = sizes[i], e = er, w = w, phi = phi), 
                        cab.rel(dist = dist, ab = b, n = sizes[i], e = er, w = w, phi = phi), 
                        set.plus.rel(dist = dist, a = a, b = b, n = sizes[i], e = er, w = w, phi = phi), 
                        set.minus.rel(dist = dist, a = a, b = b, n = sizes[i], e = er, w = w, phi = phi)
                        )
        }
      }
    }
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nSample size:\n")
  print(matrix(c(sizes, 1 - rho), length(sizes), 2, 
               dimnames = list(rep("  ", length(sizes)), c(" n", " 1-rho"))))
  cat("\n")
} 
