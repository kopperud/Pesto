
#' Title
#'
#' @param lambda
#' @param mu
#' @param eta
#' @param tstart
#' @param tend
#' @param E0
#' @param D0
#' @param ntimes
#'
#' @return
#' @export
#'
#' @examples
backwards <- function(lambda, mu, eta, tstart, tend, E0, D0, ntimes){
  k <- length(lambda)
  parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta, "k" = k)
  yini <- c("E" = E0, "D" = D0)

  times <- seq(tstart, tend, length.out = ntimes)
  out <- deSolve::rk4(y = yini, times = times, func = backwardsED, parms = parameters)

  res <- as_tibble(apply(out, 2, function(l) l, simplify = FALSE))
  return(res)
}

backwardsED <- function(t, y, parms) {
  E <- parms[["mu"]] - (parms[["lambda"]] + parms[["mu"]] + parms[["eta"]]) * y[1:parms[["k"]]] + parms[["lambda"]] * y[1:parms[["k"]]]^2 + (parms[["eta"]]/(parms[["k"]]-1)) * (sum(y[1:parms[["k"]]]) - y[1:parms[["k"]]])
  D <- - (parms[["lambda"]] + parms[["mu"]] + parms[["eta"]]) * y[(parms[["k"]]+1):(2*parms[["k"]])] + 2 * parms[["lambda"]] * y[(parms[["k"]]+1):(2*parms[["k"]])] * y[1:parms[["k"]]] + (parms[["eta"]] / (parms[["k"]] - 1)) * (sum(y[(parms[["k"]]+1):(2*parms[["k"]])]) - y[(parms[["k"]]+1):(2*parms[["k"]])])

  return(list(c("E" = E, "D" = D)))
}

