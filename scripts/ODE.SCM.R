
stochastic.character.mapping <- function(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, METHOD) {

    k <- length(lambda)
    num_branches <- length(branch_lengths)


    D_ends <- array(0,k*num_branches)
    E_ends <- array(0,k*num_branches)
    dim(D_ends) <- c(num_branches,k)
    dim(E_ends) <- c(num_branches,k)

    for (i in 1:num_branches) {

        ## calculate along the branch
        bl <- branch_lengths[i]

        ## get the initial values at the beginning of the branch for
        ## * D
        ## * E
        ## if this is a terminal branch, then we use the initial values for D and E,
        ## otherwise we need to use the previous E and the product of D_i * D_j
        if ( i %in% parents ) {
            ## this is an internal branch
            my_children <- which( parents == i )
            D_start <- lambda    ## remember lambda is a vector
            for (j in my_children) D_start <- D_start * D_ends[j,]

            ## initialize E
            E_start <- E_ends[my_children[1],]
        } else {
            ## this is an external branch
            D_start <- D_inits[i,]
            E_start <- array(0,k)
        }

        tmp_res <- branch.prob.backwards(lambda, mu, eta, bl, D_start, E_start, STEPS)
        D_ends[i,] <- tmp_res$D
        E_ends[i,] <- tmp_res$E

    }

    ## compute the root probabilities
    root_children <- which( !is.finite(parents) )
    root_probs <- lambda
    for (j in root_children) root_probs <- root_probs * D_ends[j,]

cat("Likelihood:\t\t",log(mean(root_probs)),"\n")
    root_probs <- root_probs / sum(root_probs)



    state_likelihoods <- array(0,k)

    observed_node <- 4
    bl <- branch_lengths[observed_node]
    E_start <- E_ends[observed_node,]

    if ( METHOD == "A" ) {
        F_start <- D_ends[3,]
        D_start <- D_ends[4,]

        tmp_res <- branch.prob.forwards(lambda, mu, eta, bl, F_start, D_start, E_start, STEPS, METHOD)
        state_likelihoods <- tmp_res$F

        my_children <- which( parents == observed_node )
        state_likelihoods <- lambda * state_likelihoods
        for (j in my_children) state_likelihoods <- state_likelihoods * D_ends[j,]

    } else if ( METHOD == "B" ) {
        F_start <- D_ends[3,] * lambda
        D_start <- array(0,k)

        tmp_res <- branch.prob.forwards(lambda, mu, eta, bl, F_start, D_start, E_start, STEPS, METHOD)
        state_likelihoods <- tmp_res$F

        my_children <- which( parents == observed_node )
        state_likelihoods <- lambda * state_likelihoods
        for (j in my_children) state_likelihoods <- state_likelihoods * D_ends[j,]

    } else if ( METHOD == "C" ) {
        D_start <- array(0,k)

        state_likelihoods <- array(0,k)

        for ( obs_state in 1:k ) {

            F_start <- array(0,k)
            F_start[obs_state] <- 1
            F_start <- F_start * lambda

            tmp_res <- branch.prob.forwards(lambda, mu, eta, bl, F_start, D_start, E_start, STEPS, METHOD)
            tmp_state_likelihoods <- tmp_res$F

            my_children <- which( parents == observed_node )
            tmp_state_likelihoods <- lambda * tmp_state_likelihoods
            for (j in my_children) tmp_state_likelihoods <- tmp_state_likelihoods * D_ends[j,]

            state_likelihoods <- state_likelihoods + tmp_state_likelihoods * D_ends[3,obs_state]

        }

    }

    ## normalize
    state_likelihoods <- state_likelihoods / sum(state_likelihoods)

    return ( list(root=root_probs ,node=state_likelihoods) )

}
