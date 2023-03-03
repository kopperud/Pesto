#' ML estimation of rates under the constant-rate Birth-Death process
#'
#' @param phy an object of class phylo
#' @param rho the taxon sampling fraction. Must be above 0.0 and at maximum 1.0
#'
#' @return
#' @export
#'
#' @examples
#' library(Pesto)
#'
#' data(primates)
#' rho <- 0.67
#'
#' res <- estimate_constant_bdp(primates, rho)
estimate_constant_bdp <- function(phy, rho = 1){
  phy <- treeprecompute(phy)

  JuliaCall::julia_library("Diversification")

  JuliaCall::julia_assign("phy", phy)
  JuliaCall::julia_assign("rho", rho)


  JuliaCall::julia_eval("data = Diversification.make_SSEdata2(phy, rho)")

  mlres <- JuliaCall::julia_eval("estimate_constant_bdp(data)")

  res <- list(
    "lambda" = mlres[[1]],
    "mu" = mlres[[2]]
  )

  return(res)
}
