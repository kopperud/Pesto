library(diversitree)
library(ape)

anc.state.prob.bisse <- function(phy, datafile, lambda, mu, eta) {

  # parameters = λ0, λ1, µ0, µ1, q01, q10
  pars = c(lambda, mu, eta, eta)
  #phy = read.nexus(treefile)
  tmp = read.csv(datafile,header=TRUE, row.names=1)
  states <- tmp[,"state"]
  names(states) <- rownames(tmp)

  # calculate likelihood
  lik = make.bisse(phy, states, control = list("backend", "deSolve"), strict=FALSE)
  rate = pars[1]
  num_taxa = length(states)
  #lnl <- lik(pars,root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=TRUE, intermediates=TRUE)
  lnl <- lik(pars,root=ROOT.FLAT, condition.surv=TRUE, intermediates=TRUE)
  
  D <- attr(lnl, "vals")[3:4]
  E <- attr(lnl, "vals")[1:2]
  lq <- attr(lnl, "intermediates")$lq
  
  freqs <- c(0.5, 0.5)
  nonextinct <- (1 - E)^2
  prob <- sum(freqs * D / (nonextinct * lambda)) 
  
  logL <- log(mean(D / (lambda * (1 - E)^2))) + sum(lq)
  cat("diversitree lnl =", logL, "\n")

  # now infer marginal ancestral state reconstructions
  anc_states = asr.marginal(lik, pars)

  anc_states <- t(anc_states)
  return(anc_states)
}
