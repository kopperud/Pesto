forwards <- function(lambda, mu, eta, tstart, tend, E0, D0, F0, ntimes = 100){
  k <- length(lambda)
  parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta, "k" = k)
  yini <- c("E" = E0, "D" = D0, "F" = F0)

  times <- seq(tstart, tend, length.out = ntimes)
  out <- deSolve::rk4(y = yini, times = times, func = forwardsEDF, parms = parameters)

  res <- as_tibble(apply(out, 2, function(l) l, simplify = FALSE))
  return(res)
}

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
