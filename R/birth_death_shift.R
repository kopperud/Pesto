#' Branch-specific speciation and extinction rate estimation
#'
#' @param phy
#' @param lambda
#' @param mu
#' @param eta
#' @param rho
#' @param verbose
#'
#' @return
#'
#' @export
#' @examples
#' library(Pesto)
#' library(ape)
#' library(ggtree)
#' library(tidytree)
#' library(ggplot2)
#' library(dplyr)
#'
#' data(primates)
#'
#' lambda <- c(0.1, 0.2, 0.3, 0.4, 0.20)
#' mu <- c(0.05, 0.15, 0.05, 0.15, 0.25)
#' rho <- 0.67
#' eta <- 0.008
#'
#' td <- birth_death_shift(primates, lambda, mu, eta, rho)
#'
#' th <- max(node.depth.edgelength(primates))
#' p1 <- ggtree(td, aes(color = `Speciation rate`)) +
#'              geom_tiplab(size = 3) +
#'              theme(legend.position = c(0.15, 0.8)) +
#'              xlim(c(0.0, th + 10))
#'
#' plot(p1)
birth_death_shift <- function(phy, lambda, mu, eta, rho, verbose = FALSE){
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

  res <- JuliaCall::julia_eval("birth_death_shift(model, data; verbose = verbose)")

  lambda_average <- res$lambda
  mu_average <- res$mu

  th <- max(ape::node.depth.edgelength(phy))
  df1 <- tibble::tibble("node" = 1:max(phy$edge),
                "Speciation rate" = lambda_average,
                "Extinction rate" = mu_average)
  x <- tibble::as_tibble(phy)

  df <- merge(x, df1, by = "node")
  td <- tidytree::as.treedata(df)

  return(td)
}
