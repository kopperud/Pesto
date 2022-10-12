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
  lik = make.bisse(phy, states, strict=FALSE)
  rate = pars[1]
  num_taxa = length(states)
  lnl = lik(pars,root=ROOT.GIVEN, root.p=c(0.5,0.5), condition.surv=FALSE, intermediates=TRUE)
  cat("diversitree lnl =", lnl, "\n")

  # now infer marginal ancestral state reconstructions
  anc_states = asr.marginal(lik, pars)

  anc_states <- t(anc_states)
  return(anc_states)
}
