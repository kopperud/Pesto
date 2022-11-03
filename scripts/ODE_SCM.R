# source("scripts/ODE.branch.backwards.R")
# source("scripts/ODE.branch.forwards.R")
# source("scripts/ODE.node.R")
# source("scripts/ODE.SCM.R")
# source("scripts/ODE.BiSSE.R")
# source("scripts/ODE.branch.backwards.rk4.R")

#library(deSolve)

library(BDS)

STEPS <- 10000

lambda <- c(2,1.0)
mu     <- c(0.5,0.1)
eta    <- 0.1
k      <- length(lambda)

setwd("~/projects/BDS_deterministic_map/")

if (F){
  treefile <- "data/test.tre"
  datafile <- "data/test.csv"
  phy <- read.tree(treefile)
}else{
  treefile <- "data/bears.tre"
  datafile <- "data/bears.csv"
  phy <- read.nexus(treefile)
}

 # treefile <- "data/fourtaxon.tre"
 # datafile <- "data/fourtaxon.csv"
 # phy <- read.tree(treefile)


D_inits <- matrix(0, nrow = nrow(phy$edge), ncol = 2)
df <- read.table(datafile, header = TRUE, sep = ",")
df <- df[match(phy$tip.label, df$Taxon), ]
for (i in 1:length(phy$tip.label)){
  D_inits[phy$edge[,2] == i, df$state[i]+1] <- 1
}

cat("\n")
cat("Testing marginal node state probability estimation\n")
cat("==================================================\n")
cat("\n")
cat("\n")


# results_ace_rb <- read.table(file="output/anc_states_BiSSE.log",sep="\t",head=TRUE)
# cat("Ancestral state algorithm (RevBayes)\n")
# cat("State probabilities at root:\t\t",mean(results_ace_rb[,"end_5"] == 0),mean(results_ace_rb[,"end_5"] == 1),"\n")
# cat("State probabilities at node:\t\t",mean(results_ace_rb[,"end_4"] == 0),mean(results_ace_rb[,"end_4"] == 1),"\n")



cat("\n")
cat("Ancestral state algorithm (BiSSE)\n")
results_ace_bisse <- anc.state.prob.bisse(phy, datafile, lambda, mu, eta)
print(results_ace_bisse)


cat("\n")
cat("Stochastic character mapping (B)\n")
results_scm <- traversal(phy, D_inits, lambda, mu, eta, STEPS)
print(results_scm$ASP)
