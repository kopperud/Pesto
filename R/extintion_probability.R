#' Compute the extinction probability E(t)
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

  JuliaCall::julia_eval("data = Pesto.make_SSEdata2(phy, rho);")
  JuliaCall::julia_eval("E = extinction_probability(model, data);", need_return = "Julia")

  u <- JuliaCall::julia_eval("hcat(E.u...)")
  age <- JuliaCall::julia_eval("E.t")

  return(list(
    u, age
  ))
}
