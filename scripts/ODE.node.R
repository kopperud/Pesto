
ancestral.state.probs <- function(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS) {

k <- length(lambda)
num_branches <- length(branch_lengths)

observed_node <- 4
state_likelihoods <- array(0,k)

for ( obs_state in 1:k ) {

    D_ends <- array(0,k*num_branches)
    E_ends <- array(0,k*num_branches)
    dim(D_ends) <- c(num_branches,k)
    dim(E_ends) <- c(num_branches,k)

    for (i in 1:num_branches) {

        ## calculate along the branch
        bl <- branch_lengths[i]
        if ( i %in% parents ) {
            ## this is an internal branch
            my_children <- which( parents == i )
            D_start <- lambda
            for (j in my_children) D_start <- D_start * D_ends[j,]
            E_start <- E_ends[my_children[1],]
        } else {
            ## this is an external branch
            D_start <- D_inits[i,]
            E_start <- array(0,k)
        }
        
        if ( i == observed_node ) {
            mask <- array(0,k)
            mask[obs_state] <- 1
            D_start <- D_start * mask
        }

        tmp_res <- branch.prob.backwards(lambda, mu, eta, bl, D_start, E_start, STEPS)
        D_ends[i,] <- tmp_res$D
        E_ends[i,] <- tmp_res$E
        
        
    }

    ## compute the root probabilities
    root_children <- which( !is.finite(parents) )
    root_probs <- array(1.0,k)
    root_probs <- lambda
    for (j in root_children) root_probs <- root_probs * D_ends[j,]
    
    state_likelihoods[obs_state] <- mean(root_probs)

}

## normalize
state_likelihoods <- state_likelihoods / sum(state_likelihoods)

return (state_likelihoods)

}
