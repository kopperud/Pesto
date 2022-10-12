# library(deSolve)
# 
# backwards <- function(t, state, parameters) {
#   E <- parameters[["E"]]
#   D <- parameters[["D"]]
#   with(as.list(c(state, parameters)), {
#     for (j in 1:k){
#       dE[j] <-  - (lambda[j]+mu[j]+eta)*E[j] + lambda[j] * E[j]^2
#       dD[j] <-  - (lambda[j]+mu[j]+eta)*D[j] + 2*lambda[j]*D[j] * E[j]
#       
#       for (m in 1:k){
#         if (j != m) {
#           dE[j] <- dE[j] + eta/(k-1) * E[m]
#           dD[j] <- dD[j] + eta/(k-1) * D[m]
#         }
#       }
#     }
#     list(c(dE, dD))
#   })
# }
# 
# parameters <- c(lambda = lambda, mu = mu, eta = eta)
# state      <- c(E = c(0,0), D = c(1,0))
# times      <- seq(0, 5, by = 0.01)
# 
# 
# out <- ode(y = state, times = times, func = backwards, parms = parameters)


branch.prob.backwards <- function(lambda, mu, eta, TIME, D.init, E.init, STEPS) {

  k     <- length(lambda)
  dt    <- TIME/STEPS

  prob_extinction_backward   <- array(0,k*(STEPS+1))
  prob_observed_backward     <- array(0,k*(STEPS+1))

  dim(prob_extinction_backward)  <- c(k,STEPS+1)
  dim(prob_observed_backward)    <- c(k,STEPS+1)


  prob_extinction_backward[,1]   <- E.init
  prob_observed_backward[,1]     <- D.init

  ## iterate over all time slices
  for (i in 1:STEPS) {

    ## iterate over all rate categories
    for (j in 1:k) {
        new_e <- mu[j] - (lambda[j]+mu[j]+eta)*prob_extinction_backward[j,i] + lambda[j]*prob_extinction_backward[j,i]^2
        new_d <- - (lambda[j]+mu[j]+eta)*prob_observed_backward[j,i] + 2*lambda[j]*prob_observed_backward[j,i]*prob_extinction_backward[j,i]
        for (m in 1:k) {
            if ( m != j ) {
                new_e <- new_e + eta/(k-1) * prob_extinction_backward[m,i]
                new_d <- new_d + eta/(k-1) * prob_observed_backward[m,i]
            }
        }
        prob_extinction_backward[j,i+1] <- prob_extinction_backward[j,i] + new_e*dt
        prob_observed_backward[j,i+1]   <- prob_observed_backward[j,i]   + new_d*dt
    }
  }


  rv = list(D=prob_observed_backward[,STEPS+1], E=prob_extinction_backward[,STEPS+1])

  return (rv)

}
