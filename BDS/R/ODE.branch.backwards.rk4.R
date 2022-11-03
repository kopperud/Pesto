
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
backwards <- function(lambda, mu, eta, tstart, tend, E0, D0, ntimes = 100){
  backwardsED <- function(t, y, parms) {
    E <- parms[["mu"]] - (parms[["lambda"]] + parms[["mu"]] + parms[["eta"]]) * y[1:parms[["k"]]] + parms[["lambda"]] * y[1:parms[["k"]]]^2 + (parms[["eta"]]/(parms[["k"]]-1)) * (sum(y[1:parms[["k"]]]) - y[1:parms[["k"]]])
    D <- - (parms[["lambda"]] + parms[["mu"]] + parms[["eta"]]) * y[(parms[["k"]]+1):(2*parms[["k"]])] + 2 * parms[["lambda"]] * y[(parms[["k"]]+1):(2*parms[["k"]])] * y[1:parms[["k"]]] + (parms[["eta"]] / (parms[["k"]] - 1)) * (sum(y[(parms[["k"]]+1):(2*parms[["k"]])]) - y[(parms[["k"]]+1):(2*parms[["k"]])])

    return(list(c("E" = E, "D" = D)))
  }
  k <- length(lambda)
  parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta, "k" = k)
  yini <- c("E" = E0, "D" = D0)

  times <- seq(tstart, tend, length.out = ntimes)
  out <- deSolve::rk4(y = yini, times = times, func = backwardsED, parms = parameters)
  return(as.data.frame(out))
}

# branch.prob.backwards.rk4 <- function(lambda, mu, eta, TIME, D.init, E.init){
#   backwards <- function(t, y, parms) {
#     E1 <-  parms[[2]][1] - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[1] + parms[[1]][1] * y[1]^2 + parms[[3]] * y[2]
#     E2 <-  parms[[2]][2] - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[2] + parms[[1]][2] * y[2]^2 + parms[[3]] * y[1]
#
#     D1 <- - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[3] + 2*parms[[1]][1]*y[3] * y[1] + parms[[3]] * y[4]
#     D2 <- - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[4] + 2*parms[[1]][2]*y[4] * y[2] + parms[[3]] * y[3]
#
#     return(list(c(E1, E2, D1, D2)))
#   }
#
#   parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta)
#   yini <- c(E1 = E.init[1],
#             E2 = E.init[2],
#             D1 = D.init[1],
#             D2 = D.init[2])
#   times <- seq(0, TIME, by = 0.001)
#   #times <- c(0, TIME)
#
#   out <- deSolve::rk4(y = yini, times = times, func = backwards, parms = parameters)
#   return(as.data.frame(out))
# }
#
# extinction.rk4 <- function(lambda, mu, eta, tstart, tend, u0, ntimes = 100){
#   backwards <- function(t, y, parms) {
#     E <- parms[[2]] - (parms[[1]] + parms[[2]] + parms[[3]]) * y + parms[[1]] * y^2 + (parms[[3]]/(parms[["k"]]-1)) * (sum(y) - y)
#     return(list(E))
#   }
#   k <- length(lambda)
#   parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta, "k" = k)
#   yini <- c("E" = u0)
#
#   times <- seq(tstart, tend, length.out = ntimes)
#   out <- deSolve::rk4(y = yini, times = times, func = backwards, parms = parameters)
#   return(as.data.frame(out))
# }

