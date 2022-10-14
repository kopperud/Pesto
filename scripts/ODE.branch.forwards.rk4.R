branch.prob.forwards.rk4 <- function(lambda, mu, eta, TIME, F.init, D.init, E.init){
  forwards <- function(t, y, parms) {
    E1 <-  parms[[2]][1] - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[1] + parms[[1]][1] * y[1]^2 + parms[[3]] * y[2]
    E2 <-  parms[[2]][2] - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[2] + parms[[1]][2] * y[2]^2 + parms[[3]] * y[1]
    
    D1 <- - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[3] + 2*parms[[1]][1]*y[3] * y[1] + parms[[3]] * y[4]
    D2 <- - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[4] + 2*parms[[1]][2]*y[4] * y[2] + parms[[3]] * y[3]
    
    F1 <- (-1) * - (parms[[1]][1]+parms[[2]][1]+parms[[3]])*y[5] + 2*parms[[1]][1]*y[5] * y[1] + parms[[3]] * y[6] 
    F2 <- (-1) * - (parms[[1]][2]+parms[[2]][2]+parms[[3]])*y[6] + 2*parms[[1]][2]*y[6] * y[2] + parms[[3]] * y[5]
    
    return(list(c(E1, E2, D1, D2, F1, F2)))
  }
  
  parameters <- list("lambda" = lambda, "mu" = mu, "eta" = eta)
  yini <- c(E1 = E.init[1],
            E2 = E.init[2],
            D1 = D.init[1],
            D2 = D.init[2],
            F1 = F.init[1],
            F2 = F_init[2])
  times <- seq(0, TIME, by = 0.001)
  #times <- c(0, TIME)
  
  out <- deSolve::rk4(y = yini, times = times, func = forwards, parms = parameters) 
  return(as.data.frame(out))
}