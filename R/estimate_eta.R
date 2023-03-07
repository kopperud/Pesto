#' Estimate the shift rate conditional on lambda and mu
#'
#' @param lambda a vector of speciation rates
#' @param mu a vector of extinction rates
#' @param phy an object of class phylo
#' @param lower lower limit for the numerical optimization of eta
#' @param upper upper limit for the numerical optimization of eta
#' @param eta0 starting value for eta
#'
#' @description finds the maximum-likelihood of eta conditional on the vectors lambda and mu, using gradient descent. This will typically evaluate the likelihood function around 20-60 times, and also the derivative of the likelihood function with respect to eta, so it is not instantaneous. Check if the result is stuck at the lower or upper bounds.
#'
#' @return
#' @export
#'
#' @examples
#' library(Pesto)
#' data(primates)
#'
#' lambda <- c(0.1, 0.2, 0.3, 0.4)
#' mu <- c(0.05, 0.15, 0.05, 0.15)
#' rho <- 0.67
#'
#' etaml <- optimize_eta(lambda, mu, primates, rho)
#' print(etaml)
optimize_eta <- function(lambda, mu, phy, rho, lower = 0.00001, upper = 10.0, eta0 = NULL){
  phy <- treeprecompute(phy)

  # some reasonable starting value if not specified
  if(is.null(eta0)){
    treelength <- sum(phy$edge.length)
    eta0 <- 1 / treelength
  }

  JuliaCall::julia_library("Pesto")

  JuliaCall::julia_assign("phy", phy)
  JuliaCall::julia_assign("lambda", lambda)
  JuliaCall::julia_assign("mu", mu)
  JuliaCall::julia_assign("rho", rho)
  JuliaCall::julia_assign("lower", lower)
  JuliaCall::julia_assign("upper", upper)
  JuliaCall::julia_assign("eta0", eta0)

  JuliaCall::julia_eval("data = Pesto.make_SSEdata2(phy, rho)")

  eta <- JuliaCall::julia_eval("optimize_eta(lambda, mu, data; lower = lower, upper = upper, xinit = eta0)")
  return(eta)
}
