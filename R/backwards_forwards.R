#' Title
#'
#' @param phy
#' @param lambda
#' @param mu
#' @param rho
#'
#' @return
#' @export
#'
#' @examples
extinction_probability <- function(phy, lambda, mu, rho){
  phy <- treeprecompute(phy)

  JuliaCall::julia_assign("phy", phy)
  JuliaCall::julia_assign("lambda", lambda)
  JuliaCall::julia_assign("mu", mu)
  JuliaCall::julia_assign("rho", rho)

  JuliaCall::julia_eval("data = Diversification.make_SSEdata2(phy, rho);")
  JuliaCall::julia_eval("E = extinction_probability(model, data);", need_return = "Julia")

  u <- JuliaCall::julia_eval("hcat(E.u...)")
  age <- JuliaCall::julia_eval("E.t")

  return(list(
    u, age
  ))
}

backwards <- function(phy, lambda, mu, rho, E){
  return(1)
  #D_ends, Ds, sf = postorder(model, data, E; verbose = verbose, alg = alg)
}

forwards <- function(phy, lambda, mu, rho, E, Ds){
  return(1)
}

# Ds, Fs = Diversification.backwards_forwards_pass(model, data);
# Ss = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)
# E = extinction_probability(model, data)
