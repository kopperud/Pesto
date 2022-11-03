#' Title
#'
#' @param phy
#' @param D_inits
#' @param lambda
#' @param mu
#' @param eta
#' @param ntimes
#'
#' @return
#' @export
#'
#' @examples
traversal <- function(phy, D_inits, lambda, mu, eta, ntimes = 100, printlogL = FALSE) {
    po <- postorder(phy)
    k <- length(lambda)
    num_branches <- length(phy$edge.length)

    forward_results <- list()

    D_ends <- matrix(0, nrow = num_branches, ncol = k)
    E_ends <- matrix(0, nrow = num_branches, ncol = k)

    sf <- matrix(0, nrow = num_branches, ncol = 1)

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

      tmp_res <- backwards(lambda, mu, eta, tstart = 0, tend = bl, E_start, D_start, ntimes = 100)
      D_end <- tmp_res %>%
        dplyr::select(dplyr::starts_with("D")) %>%
        tail(n = 1) %>%
        unlist()

      E_ends[i,] <- tmp_res %>%
        dplyr::select(dplyr::starts_with("E")) %>%
        tail(n = 1) %>%
        unlist()

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

    if (printlogL){
      # probability of at least two lineages not going extinct before the present
      nonextinct <- (1 - E_ends[root_children[1],]) ^2
      logL <- log(mean(root_probs / (nonextinct * lambda))) + sum(sf)

      cat("Likelihood:\t\t", logL,"\n")
    }
    root_probs <- root_probs / sum(root_probs)

    #############################################
    ##
    ##      Preorder pass, compute `F(t)`
    ##
    #############################################
    root_node <- length(phy$tip.label) + 1
    left_root_edge <- which(phy$edge[,1] == root_node)[1]

    E_start <- E_ends[left_root_edge,]

    F_ends <- matrix(0, nrow = nrow(phy$edge), ncol = k)

    ## Preorder traversal
    for (i in rev(po)){
      parent <- phy$edge[i,1]
      child <- phy$edge[i,2]
      bl <- phy$edge.length[i]

      ## if root
      if(parent == length(phy$tip.label) + 1){
        root_children <- which(phy$edge[,1] == length(phy$tip.label) +1)
        other_child <- root_children[root_children != i]
        F_start <- D_ends[other_child,] * lambda
      }else{
        parent_edge <- which(phy$edge[,2] == parent)
        children <- which(phy$edge[,1] == parent)
        other_child <- children[children != i]

        F_start <- F_ends[parent_edge,] * lambda * D_ends[other_child,]

        ## Normalize
        F_start <- F_start / sum(F_start)
      }

      D_start <- D_ends[i,]
      E_start <- E_ends[i,]


      tmp_res <- forwards(lambda, mu, eta, tstart = 0, tend = -bl, E_start, D_start, F_start, ntimes = 100)
      forward_results[[i]] <- forwards

      F_ends[i,] <- tmp_res %>%
        dplyr::select(dplyr::starts_with("F")) %>%
        tail(n = 1) %>%
        unlist()

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
    }

#       x <- seq(0, bl, length.out = STEPS+1)
#       df1 <- dplyr::tibble("time" = x,
#                    "F" = tmp_res$F[1,],
#                    "state" = "1")
#       df2 <- dplyr::tibble("time" = x,
#                    "F" = tmp_res$F[2,],
#                    "state" = "2")
#       df <- dplyr::bind_rows(df1, df2)
#       p <- ggplot2::ggplot(df, ggplot2::aes(x = time, y = F, linetype = state)) +
#                   ggplot2::geom_line()
#       ggplot2::ggsave("figures/tmp.pdf", p)


    res <- list(
      "ASP" = ASP,
      "forwards" = forward_results
    )
    return(res)
}
