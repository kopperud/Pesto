
#stochastic.character.mapping <- function(branch_lengths, parents, D_inits, lambda, mu, eta, STEPS, METHOD) {
stochastic.character.mapping <- function(phy, D_inits, lambda, mu, eta, STEPS, METHOD) {

    po <- postorder(phy)
    k <- length(lambda)
    num_branches <- length(phy$edge.length)


    D_ends <- array(0,k*num_branches)
    E_ends <- array(0,k*num_branches)
    dim(D_ends) <- c(num_branches,k)
    dim(E_ends) <- c(num_branches,k)
    
    sf <- matrix(0, nrow = length(phy$edge), ncol = 1)

    for (i in po){
      parent <- phy$edge[i,1]
      child <- phy$edge[i,2]
      

        ## calculate along the branch
        bl <- phy$edge.length[i]

        ## get the initial values at the beginning of the branch for
        ## * D
        ## * E
        ## if this is a terminal branch, then we use the initial values for D and E,
        ## otherwise we need to use the previous E and the product of D_i * D_j
        if ( child > length(phy$tip.label) ) {
            ## this is an internal branch
            my_children <- which(phy$edge[,1] == child)
            D_start <- lambda
            for (j in my_children){
              D_start <- D_start * D_ends[j,]
            } 

            ## initialize E
            E_start <- E_ends[my_children[1],]
        } else {
            ## this is an external branch
            D_start <- D_inits[i,]
            E_start <- array(0,k)
        }

        
        if(FALSE){
          tmp_res <- branch.prob.backwards(lambda, mu, eta, bl, D_start, E_start, STEPS)
          D_end <- tmp_res$D
          E_ends[i,] <- tmp_res$E
        }else{
          tmp_res <- branch.prob.backwards.rk4(lambda, mu, eta, bl, D_start, E_start)
          D_end <- c(tail(tmp_res$D1, n = 1), tail(tmp_res$D2, n = 1))

          E_ends[i,] <- c(tail(tmp_res$E1, n = 1), tail(tmp_res$E2, n = 1))
        }
        sf[i,] <- log(sum(D_end)) ## add the scaling factor
        D_end <- D_end / (sum(D_end))
        D_ends[i,] <- D_end
      }

    
    ## compute the root probabilities
    root_node <- length(phy$tip.label) + 1
    root_children <- which(phy$edge[,1] == root_node)
    root_probs <- lambda
    for (j in root_children){
      root_probs <- root_probs * D_ends[j,]
    }

    #cat("Likelihood:\t\t",log(mean(root_probs)) + sum(sf),"\n")
    # probability of at least two lineages not going extinct before the present
    nonextinct <- (1 - E_ends[root_children[1],]) ^2
    logL <- log(mean(root_probs / (nonextinct * lambda))) + sum(sf)
    
    cat("Likelihood:\t\t", logL,"\n")
    root_probs <- root_probs / sum(root_probs)

    # n <- length(phy$tip.label); (n-1) * log(2) - sum(log(1:n)) - sum(log(1:(n-1)))

    state_likelihoods <- array(0,k)

    #############################################
    ##
    ##      Preorder pass, compute `F(t)`
    ## 
    #############################################
    root_node <- length(phy$tip.label) + 1
    left_root_edge <- which(phy$edge[,1] == root_node)[1]
    
    E_start <- E_ends[left_root_edge,]

    if ( METHOD == "A" ) {
        F_start <- D_ends[3,]
        D_start <- D_ends[4,]

        tmp_res <- branch.prob.forwards(lambda, mu, eta, bl, F_start, D_start, E_start, STEPS, METHOD)
        state_likelihoods <- tmp_res$Froot

        my_children <- which( parents == observed_node )
        state_likelihoods <- lambda * state_likelihoods
        for (j in my_children) state_likelihoods <- state_likelihoods * D_ends[j,]

    } else if ( METHOD == "B" ) {
        F_ends <- matrix(0, nrow = nrow(phy$edge), ncol = k)
        
        for (i in rev(po)){
          parent <- phy$edge[i,1]
          child <- phy$edge[i,2]
          bl <- phy$edge.length[i]
          
          ## if root
          if(parent == length(phy$tip.label) + 1){
            root_children <- which(phy$edge[,1] == length(phy$tip.label) +1)
            other_child <- root_children[root_children != i]
            #F_start <- D_ends[root_children[1],] * D_ends[root_children[2],] * lambda
            F_start <- D_ends[other_child,] * lambda
          }else{
            parent_edge <- which(phy$edge[,2] == parent)
            children <- which(phy$edge[,1] == parent)
            other_child <- children[children != i]
            
            F_start <- F_ends[parent_edge,] * lambda * D_ends[other_child,]
          }
          
          D_start <- D_ends[i,] 
          E_start <- E_ends[i,]
        
          tmp_res <- branch.prob.forwards(lambda, mu, eta, bl, F_start, D_start, E_start, STEPS, METHOD)
          F_ends[i,] <- tmp_res$Froot
        }

        ASP <- matrix(0, nrow = length(phy$tip.label) - 1, ncol = k)
        
        for (node in (length(phy$tip.label) + 1): max(phy$edge)){
          children <- which(phy$edge[,1] == node)
          
          if (node == length(phy$tip.label) + 1){
            ASP[node - length(phy$tip.label), ] <- D_ends[children[1],] * D_ends[children[2],] * lambda
          }else{
            edge_idx <- which(phy$edge[,2] == node)
            ASP[node - length(phy$tip.label), ] <- D_ends[children[1],] * D_ends[children[2],] * lambda * F_ends[edge_idx,]
          }
          ASP[node - length(phy$tip.label),] <- ASP[node - length(phy$tip.label),] / sum(ASP[node - length(phy$tip.label),])
        }
        
        
    } else if ( METHOD == "C" ) {
        D_start <- array(0,k)

        state_likelihoods <- array(0,k)

        for ( obs_state in 1:k ) {

            F_start <- array(0,k)
            F_start[obs_state] <- 1
            F_start <- F_start * lambda

            print("F_start:")
            print(F_start)

            tmp_res <- branch.prob.forwards(lambda, mu, eta, bl, F_start, D_start, E_start, STEPS, METHOD)
            tmp_state_likelihoods <- tmp_res$Froot

            my_children <- which( parents == observed_node )
            tmp_state_likelihoods <- lambda * tmp_state_likelihoods
            for (j in my_children) tmp_state_likelihoods <- tmp_state_likelihoods * D_ends[j,]

            state_likelihoods <- state_likelihoods + tmp_state_likelihoods * D_ends[3,obs_state]
        }
        
      
      

      x <- seq(0, bl, length.out = STEPS+1)
      df1 <- dplyr::tibble("time" = x,
                   "F" = tmp_res$F[1,],
                   "state" = "1")
      df2 <- dplyr::tibble("time" = x,
                   "F" = tmp_res$F[2,],
                   "state" = "2")
      df <- dplyr::bind_rows(df1, df2)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = F, linetype = state)) +
                  ggplot2::geom_line()
      ggplot2::ggsave("figures/tmp.pdf", p)


    }

    ## normalize
    # state_likelihoods <- state_likelihoods / sum(state_likelihoods)
    

    res <- list(
      "ASP" = ASP
    )
    # return ( list(root=root_probs ,node=state_likelihoods) )

}
