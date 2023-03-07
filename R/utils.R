#' Precompute tree statistics
#'
#' @param phy
#'
#' @description Add some tree statistics to the tree object, compute the branch traversal order (postorder).
#'
#' @return
#' @export
#'
#' @examples
#' library(Pesto)
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

#' lognormal quantiles
#'
#' @param meanlog the mean of the normal distribution on a log scale
#' @param sdlog the standard deviation of the normal distribution on a log scale. sdlog = 0.587 represents a distribution whose 2.5-97.5% quantile spans one order of magnitude
#' @param k the number of discrete quantiles to be returned
#'
#' @return
#' @export
#'
#' @examples
#' # Eight points from the discretized log-normal distribution
#' # with median = 0.2 and standard deviation H such that the
#' # 2.5-97.5% quantile spans one order of magnitude.
#'
#' H <- 0.587
#' q <- lognorm_quantiles(log(0.2), sdlog = H, n = 8)
lognorm_quantiles <- function(meanlog = log(0.2), sdlog = 0.587, n = 6){
  ns <- 1:n
  p <- ((ns-0.5)/n)
  q <- qlnorm(p, meanlog = meanlog, sdlog = sdlog)
  return(q)
}

#' All pairwise combinations among x and y
#'
#' @param ys
#' @param xs
#'
#' @return All the pairwise combinations
#' @export
#'
#' @examples
#' n <- 6
#' H <- 0.578
#'
#' lambda_quantiles <- lognorm_quantiles(meanlog = log(0.28), sdlog = H, n = n)
#' mu_quantiles <- lognorm_quantiles(meanlog = log(0.20), sdlog = H, n = n)
#'
#' allpairwise(lambda_quantiles, mu_quantiles)
allpairwise <- function(xs, ys){
  nx <- length(xs)
  ny <- length(ys)

  k <- ny * nx

  lambda <- numeric(k)
  mu <- numeric(k)

  i <- 1.0
  for (y in ys){
    for (x in xs){
      lambda[i] <- x
      mu[i] <- y
      i <- i + 1
    }
  }
  res <- list(
    lambda,
    mu
  )
  return(res)
}



