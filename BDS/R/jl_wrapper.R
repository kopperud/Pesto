
#' @inherit birth_death_shift
#'
#' @export
#' @examples
#' library(BDS)
#' library(ape)
#' library(ggtree)
#' library(tidytree)
#' library(ggplot2)
#' library(dplyr)
#'
#' data(bears)
#'
#' lambda <- c(0.1, 0.2)
#' mu <- c(0.05, 0.15)
#' rho <- 1.0
#' eta <- 0.05
#'
#' res <- birth_death_shift3(bears, lambda, mu, eta, rho)
#' lambda_average <- res$lambda
#'
#' th <- max(node.depth.edgelength(bears))
#' df1 <- tibble("node" = 1:max(bears$edge),
#'               "Speciation rate" = lambda_average)
#' x <- as_tibble(bears)
#'
#' bearsdf <- merge(x, df1, by = "node")
#' td_bears <- as.treedata(bearsdf)
#'
#' p1a <- ggtree(td_bears, aes(color = `Speciation rate`)) +
#'   geom_tiplab(size = 8) +
#'   theme(legend.position = c(0.2, 0.8)) +
#'   xlim(c(0.0, th + 10))
birth_death_shift3 <- function(phy, lambda, mu, eta, rho, verbose = FALSE){
  phy <- treeprecompute(phy)

  JuliaCall::julia_library("Diversification")

  JuliaCall::julia_assign("phy", phy)
  JuliaCall::julia_assign("lambda", lambda)
  JuliaCall::julia_assign("mu", mu)
  JuliaCall::julia_assign("eta", eta)
  JuliaCall::julia_assign("rho", rho)
  JuliaCall::julia_assign("verbose", verbose)

  JuliaCall::julia_eval("data = Diversification.make_SSEdata2(phy, rho)")
  JuliaCall::julia_eval("model = Diversification.SSEconstant(lambda, mu, eta)")

  avg_node_rates <- JuliaCall::julia_eval("birth_death_shift(model, data; verbose = verbose)")
  return(avg_node_rates)
}

#' Title
#'
#' @param phy
#'
#' @return
#' @export
#'
#' @examples
#' library(BDS)
#'
#' data(bears)
#'
#' bears <- treeprecompute(bears)
treeprecompute <- function(phy){
  nde <- ape::node.depth.edgelength(phy)
  node_depths <- max(nde) - nde

  phy$node_depths <- node_depths
  phy$branching_times <- ape::branching.times(phy)

  po <- ape::postorder(phy)
  phy$po <- po
  return(phy)
}
