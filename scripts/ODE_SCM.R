source("scripts/ODE.branch.backwards.R")
source("scripts/ODE.branch.forwards.R")
source("scripts/ODE.node.R")
source("scripts/ODE.SCM.R")
source("scripts/ODE.BiSSE.R")


STEPS <- 100000

lambda <- c(2,1.0)
mu     <- c(0.5,0.1)
eta    <- 0.1
k      <- length(lambda)


if (TRUE){
  treefile <- "data/test.tre"  
  datafile <- "data/test.csv"
  phy <- read.tree(treefile)
}else{
  treefile <- "data/bears.tre"
  datafile <- "data/bears.csv"
  phy <- read.nexus(treefile)  
}
 
 treefile <- "data/fourtaxon.tre"  
 datafile <- "data/fourtaxon.csv"
 phy <- read.tree(treefile)


D_inits <- matrix(0, nrow = nrow(phy$edge), ncol = 2)
df <- read.table(datafile, header = TRUE, sep = ",")
df <- df[match(phy$tip.label, df$Taxon), ]
for (i in 1:length(phy$tip.label)){
  D_inits[phy$edge[,2] == i, df$state[i]+1] <- 1
}

# branch_lengths <- c( 1.0, 1.0, 2.0, 1.0 )
# parents <- c( 4, 4, NA, NA )
# D_inits <- c( 0, 1, 0, 1, 0, 1 )
# dim(D_inits) <- c(3,k)


cat("\n")
cat("Testing marginal node state probability estimation\n")
cat("==================================================\n")
cat("\n")
cat("\n")


results_ace_rb <- read.table(file="output/anc_states_BiSSE.log",sep="\t",head=TRUE)
cat("Ancestral state algorithm (RevBayes)\n")
cat("State probabilities at root:\t\t",mean(results_ace_rb[,"end_5"] == 0),mean(results_ace_rb[,"end_5"] == 1),"\n")
cat("State probabilities at node:\t\t",mean(results_ace_rb[,"end_4"] == 0),mean(results_ace_rb[,"end_4"] == 1),"\n")



cat("\n")
cat("Ancestral state algorithm (BiSSE)\n")
results_ace_bisse <- anc.state.prob.bisse(phy, datafile, lambda, mu, eta)
#cat("State probabilities at root:\t\t",results_ace_bisse$root,"\n")
#cat("State probabilities at node:\t\t",results_ace_bisse$node,"\n")
print(results_ace_bisse)



#cat("\n")
#cat("Ancestral state algorithm\n")
#results_ace <- ancestral.state.probs(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS)
#cat("State probabilities at node:\t\t",results_ace,"\n")




#results_scm <- stochastic.character.mapping(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, "A")
#
#cat("Stochastic character mapping (A)\n")
#cat("State probabilities at root:\t\t",results_scm$root,"\n")
#cat("State probabilities at node:\t\t",results_scm$node,"\n")
#



results_scm <- stochastic.character.mapping(phy, D_inits, lambda, mu, eta, STEPS, "B")

cat("Stochastic character mapping (B)\n")
# cat("State probabilities at root:\t\t",results_scm$root,"\n")
# cat("State probabilities at node:\t\t",results_scm$node,"\n")
print(results_scm$ASP)

q()

cat("\n")
cat("Stochastic character mapping (C)\n")
results_scm <- stochastic.character.mapping(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, "C")
cat("State probabilities at root:\t\t",results_scm$root,"\n")
cat("State probabilities at node:\t\t",results_scm$node,"\n")
