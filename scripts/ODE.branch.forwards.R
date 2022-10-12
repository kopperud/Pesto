

branch.prob.forwards <- function(lambda, mu, eta, TIME, F.init, D.init, E.init, STEPS, METHOD) {

  k     <- length(lambda)
  dt    <- TIME/STEPS

  prob_extinction_forward   <- array(0,k*(STEPS+1))
  prob_observed_forward     <- array(0,k*(STEPS+1))
  prob_observed_backward    <- array(0,k*(STEPS+1))

  dim(prob_extinction_forward)  <- c(k,STEPS+1)
  dim(prob_observed_forward)    <- c(k,STEPS+1)
  dim(prob_observed_backward)   <- c(k,STEPS+1)


  prob_extinction_forward[,STEPS+1]   <- E.init
  prob_observed_forward[,STEPS+1]     <- F.init
  prob_observed_backward[,STEPS+1]    <- D.init

  ## iterate over all time slices
  for (i in STEPS:1) {

    ## iterate over all rate categories
    for (j in 1:k) {
        new_e <- mu[j] - (lambda[j]+mu[j]+eta)*prob_extinction_forward[j,i+1] + lambda[j]*prob_extinction_forward[j,i+1]^2
        new_f <- - (lambda[j]+mu[j]+eta)*prob_observed_forward[j,i+1]  + 2*lambda[j]*prob_observed_forward[j,i+1] *prob_extinction_forward[j,i+1]
        new_d <- - (lambda[j]+mu[j]+eta)*prob_observed_backward[j,i+1] + 2*lambda[j]*prob_observed_backward[j,i+1]*prob_extinction_forward[j,i+1]
        for (m in 1:k) {
            if ( m != j ) {
                new_e <- new_e + eta/(k-1) * prob_extinction_forward[m,i+1]
                new_d <- new_d + eta/(k-1) * prob_observed_backward[m,i+1]
                new_f <- new_f + eta/(k-1) * prob_observed_forward[m,i+1]
            }
        }
        prob_extinction_forward[j,i] <- prob_extinction_forward[j,i+1] - new_e*dt
        prob_observed_backward[j,i]  <- prob_observed_backward[j,i+1]  - new_d*dt
        if ( METHOD == "A" ) {
            prob_observed_forward[j,i]   <- prob_observed_forward[j,i+1]   + new_f*dt * prob_observed_backward[j,i]
        } else {
            prob_observed_forward[j,i]   <- prob_observed_forward[j,i+1]   + new_f*dt
        }
    }
  }

  rv = list("D" = prob_observed_backward, 
            "E" = prob_extinction_forward,
            "F" = prob_observed_forward,
            "Droot" = prob_observed_backward[,1], 
            "Eroot" = prob_extinction_forward[,1],
            "Froot" = prob_observed_forward[,1])

  return (rv)


}
