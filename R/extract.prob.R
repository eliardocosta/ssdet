#' Extract the probability distribution of an object of the function "exp.postDPmix".
#'
#' @param obj An object of the function "exp.postDPmix".
#'
#' @return A list with the posterior probabilities.
#'
#' @noRd
#' 
extract.prob <- function(obj) { 
  grid.val <- filter(obj$lam, rep(1/2, 2)) 
  grid.val <- as.numeric(head(grid.val, -1))
  prob <- numeric()
  for (i in 2:length(obj$cumprob)) {
    prob[i] <- obj$cumprob[i]-obj$cumprob[i-1]
  }
  prob <- prob[-1]
  prob <- prob/sum(prob)
  return(list(lam = grid.val, prob = prob))
} 
