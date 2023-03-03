#' Compute the expected number of shifts in rate category per branch
#'
#' @param phy
#' @param lambda
#' @param mu
#' @param rho
#' @param ntimeslices
#'
#' @return
#' @export
#'
#' @examples
#' library(Pesto)
#' library(ggtree)
#' library(ggplot2)
#' data(primates)
#'
#' lambda <- c(0.1, 0.2, 0.3, 0.4)
#' mu <- c(0.05, 0.15, 0.05, 0.15)
#' eta <- 0.008
#' rho <- 0.67
#'
#' td <- compute_nshifts(primates, lambda, mu, eta, rho)
#'
#' th <- max(ape::node.depth.edgelength(primates))
#' p1 <- ggtree(td, aes(color = nshifts)) +
#'              geom_tiplab(size = 3) +
#'              theme(legend.position = c(0.15, 0.8)) +
#'              xlim(c(0.0, th + 10)) +
#'              scale_colour_gradient(low = "black", high = "red")
#'
#' plot(p1)
compute_nshifts <- function(phy, lambda, mu, eta, rho, ntimeslices = 500){
  phy <- treeprecompute(phy)

  JuliaCall::julia_assign("phy", phy)
  JuliaCall::julia_assign("lambda", lambda)
  JuliaCall::julia_assign("mu", mu)
  JuliaCall::julia_assign("eta", eta)
  JuliaCall::julia_assign("rho", rho)
  JuliaCall::julia_assign("ntimeslices", ntimeslices)

  JuliaCall::julia_library("Diversification")

  JuliaCall::julia_eval("data = Diversification.make_SSEdata2(phy, rho)")
  JuliaCall::julia_eval("model = Diversification.SSEconstant(lambda, mu, eta)")

  JuliaCall::julia_eval("E = extinction_probability(model, data)", need_return="Julia")
  JuliaCall::julia_eval("Ds, Fs = Diversification.backwards_forwards_pass(model, data)", need_return="Julia")
  JuliaCall::julia_eval("Ss = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)", need_return="Julia")

  nshifts <- JuliaCall::julia_eval("nshifts = compute_nshifts(model, data, Ds, Ss; ntimeslices = ntimeslices)")

  th <- max(ape::node.depth.edgelength(phy))
  df1 <- tibble::tibble("node" = 1:max(phy$edge),
                        "nshifts" = nshifts)
  x <- tibble::as_tibble(phy)

  df <- merge(x, df1, by = "node")
  td <- tidytree::as.treedata(df)

  return(td)
}
