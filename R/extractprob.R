#' Extract the probability distribution of an object of the function "exp.postDPmix".
#'
#' @param obj An object of the function "exp.postDPmix".
#'
#' @return A probability distribution.
#'
#' @noRd
#' 
extract.prob <- function(obj) { 
  grid.val <- filter(obj$lam, rep(1/2, 2)) 
  grid.val <- as.numeric(head(grid.val, -1))
  dprob <- numeric()
  for (i in 2:length(obj$fdist)) {
    dprob[i] <- obj$fdist[i]-obj$fdist[i-1]
  }
  dprob <- dprob[-1]
  dprob <- dprob/sum(dprob)
  return(list(lam = grid.val, ddist = dprob))
} 
