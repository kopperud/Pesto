source("scripts/ODE.branch.backwards.R")
source("scripts/ODE.branch.forwards.R")
source("scripts/ODE.node.R")
source("scripts/ODE.SCM.R")
source("scripts/ODE.BiSSE.R")


STEPS <- 10000

lambda <- c(2,1.0)
mu     <- c(0.5,0.1)
eta    <- 0.1
k      <- length(lambda)

branch_lengths <- c( 1.0, 1.0, 2.0, 1.0 )
parents <- c( 4, 4, NA, NA )
D_inits <- c( 0, 1, 0, 1, 0, 1 )
dim(D_inits) <- c(3,k)


cat("\n")
cat("Testing marginal node state probability estimation\n")
cat("==================================================\n")
cat("\n")
cat("\n")


results_ace_rb <- read.table(file="output/anc_states_BiSSE.log",sep="\t",head=TRUE)
cat("Ancestral state algorithm (RevBayes)\n")
cat("State probabilities at root:\t\t",mean(results_ace_rb[,"end_5"] == 0),mean(results_ace_rb[,"end_5"] == 1),"\n")
cat("State probabilities at node:\t\t",mean(results_ace_rb[,"end_4"] == 0),mean(results_ace_rb[,"end_4"] == 1),"\n")



#results_scm_rb <- read.table(file="output/events.tsv",sep="\t",head=TRUE)
#cat("Stochastic character mapping algorithm (RevBayes)\n")



#cat("State probabilities at root:\t\t",mean(results_scm_rb[results_scm_rb[,"node_index"] == 5,"end_state"] == 0),mean(results_scm_rb[results_scm_rb[,"node_index"] == 5,"end_state"] == 1),"\n")
#cat("State probabilities at root:\t\t",mean(results_scm_rb[results_scm_rb[,"node_index"] == 5,"start_state"] == 0),mean(results_scm_rb[results_scm_rb[,"node_index"] == 5,"start_state"] == 1),"\n")
#cat("State probabilities at node:\t\t",mean(results_scm_rb[results_scm_rb[,"node_index"] == 4,"end_state"] == 0),mean(results_scm_rb[results_scm_rb[,"node_index"] == 4,"end_state"] == 1),"\n")
#cat("State probabilities at node:\t\t",mean(results_scm_rb[results_scm_rb[,"node_index"] == 4,"start_state"] == 0),mean(results_scm_rb[results_scm_rb[,"node_index"] == 4,"start_state"] == 1),"\n")

#tmp <- results_scm_rb[results_scm_rb[,"node_index"] == 4,]
#tmp2 <- tmp[with(tmp, order(transition_time, decreasing=TRUE)),]
#tmp2 <- tmp
#tmp3 <- tmp2[duplicated(tmp2[,"iteration"],fromLast=TRUE) == FALSE,]

#cat("State probabilities at node:\t\t",mean(tmp3[, "end_state"] == 0),mean(tmp3[, "end_state"] == 1),"\n")

#num_iterations = max(results_scm_rb$iteration)
#num_root_0 = 0
#num_node_0 = 0
#for (i in 0:num_iterations) {
#
#    iteration_data = results_scm_rb[results_scm_rb$iteration == i,]
#
#    # dont need to check for anagenetic changes on the root because it has no branch length
#    d = iteration_data[iteration_data$end_state == 0 & iteration_data$node_index == 5,]
#    num_root_0 = num_root_0 + nrow(d)
#
#    # for the internal nodes we must check for anagenetic change
#    d = iteration_data[iteration_data$transition_type == "anagenetic" & iteration_data$node_index == 4,]
#    if (nrow(d) > 0) {
#        # count only the last anagenetic change
#        d = head(d[with(d, order(transition_time, decreasing=!TRUE)),], 1)
#        num_node_0 = num_node_0 + (d$end_state == 0)
#    } else {
#        # no anagenetic change
#        d = iteration_data[iteration_data$end_state == 0 & iteration_data$node_index == 4,]
#        num_node_0 = num_node_0 + nrow(d)
#    }
#}
#cat("\n\nStoch map probs:\n")
#cat("State probabilities at root:\t\t",num_root_0/(num_iterations + 1), 1-(num_root_0/(num_iterations + 1)),"\n")
#cat("State probabilities at node:\t\t",num_node_0/(num_iterations + 1), 1-(num_node_0/(num_iterations + 1)),"\n")

results_ace_bisse <- anc.state.prob.bisse(lambda, mu, eta)

cat("Ancestral state algorithm (BiSSE)\n")
cat("State probabilities at root:\t\t",results_ace_bisse$root,"\n")
cat("State probabilities at node:\t\t",results_ace_bisse$node,"\n")



results_ace <- ancestral.state.probs(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS)

cat("Ancestral state algorithm\n")
cat("State probabilities at node:\t\t",results_ace,"\n")




results_scm <- stochastic.character.mapping(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, "A")

cat("Stochastic character mapping (A)\n")
cat("State probabilities at root:\t\t",results_scm$root,"\n")
cat("State probabilities at node:\t\t",results_scm$node,"\n")


results_scm <- stochastic.character.mapping(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, "B")

cat("Stochastic character mapping (B)\n")
cat("State probabilities at root:\t\t",results_scm$root,"\n")
cat("State probabilities at node:\t\t",results_scm$node,"\n")


results_scm <- stochastic.character.mapping(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, "C")

cat("Stochastic character mapping (C)\n")
cat("State probabilities at root:\t\t",results_scm$root,"\n")
cat("State probabilities at node:\t\t",results_scm$node,"\n")
