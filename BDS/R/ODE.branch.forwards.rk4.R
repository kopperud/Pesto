forwards <- function(lambda, mu, eta, tstart, tend, E0, D0, F0, ntimes = 100){
  forwardsEDF <- function(t, y, parms) {
    lambda <- parms[["lambda"]]
    mu <- parms[["mu"]]
    eta <- parms[["eta"]]
    k <- parms[["k"]]

    E <- (-1) * (mu - (lambda + mu + eta) * y[1:k] + lambda * y[1:k]^2 + (eta/(k-1)) * (sum(y[1:k]) - y[1:k]))
    D <- (-1) * (- (lambda + mu + eta) * y[(k+1):(2*k)] + 2 * lambda * y[(k+1):(2*k)] * y[1:k] + (eta / (k-1)) * (sum(y[(k+1):(2*k)]) - y[(k+1):(2*k)]))
    F1 <- (-1) * (- (lambda + mu + eta) * y[(2*k+1):(3*k)] + 2 * lambda * y[(2*k+1):(3*k)] * y[1:k] + (eta / (k-1)) * (sum(y[(2*k+1):(3*k)]) - y[(2*k+1):(3*k)]))

    return(list(c("E" = E, "D" = D, "F" = F1)))
  }

  k <- length(lambda)
  parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta, "k" = k)
  yini <- c("E" = E0, "D" = D0, "F" = F0)

  times <- seq(tstart, tend, length.out = ntimes)
  out <- deSolve::rk4(y = yini, times = times, func = forwardsEDF, parms = parameters)
  return(as.data.frame(out))
}

# branch.prob.forwards.rk4 <- function(lambda, mu, eta, TIME, F.init, D.init, E.init){
#   forwards <- function(t, y, parms) {
#     E1 <-  parms[[2]][1] - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[1] + parms[[1]][1] * y[1]^2 + parms[[3]] * y[2]
#     E2 <-  parms[[2]][2] - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[2] + parms[[1]][2] * y[2]^2 + parms[[3]] * y[1]
#
#     D1 <- - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[3] + 2*parms[[1]][1]*y[3] * y[1] + parms[[3]] * y[4]
#     D2 <- - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[4] + 2*parms[[1]][2]*y[4] * y[2] + parms[[3]] * y[3]
#
#     F1 <- (-1) * - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[5] + 2*parms[[1]][1]*y[5] * y[1] + parms[[3]] * y[6]
#     F2 <- (-1) * - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[6] + 2*parms[[1]][2]*y[6] * y[2] + parms[[3]] * y[5]
#
#     return(list(c(E1, E2, D1, D2, F1, F2)))
#   }
#
#   parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta)
#   yini <- c(E1 = E.init[1],
#             E2 = E.init[2],
#             D1 = D.init[1],
#             D2 = D.init[2],
#             F1 = F.init[1],
#             F2 = F_init[2])
#   times <- seq(0, TIME, by = 0.001)
#   #times <- c(0, TIME)
#
#   out <- deSolve::rk4(y = yini, times = times, func = forwards, parms = parameters)
#   return(as.data.frame(out))
# }


