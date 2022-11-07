#' Title
#'
#' @param phy an object of class phylo
#' @param lambda vector of speciation rates
#' @param mu vector of extinction rates
#' @param eta the (common) rate of change between state i and state j. A scalar
#' @param ntimes number of solutions provided by the ODE solver for each branch
#'
#' @return
#' @export
#'
#' @examples
birth_death_shift <- function(phy, lambda, mu, eta, ntimes = 100){
  k <- length(lambda)

  stopifnot(length(mu) == length(lambda))

  D_inits <- matrix(1, nrow = nrow(phy$edge), ncol = k)

  x <- traversal(phy, D_inits, lambda, mu, eta, ntimes = ntimes)

  n_edges <- nrow(phy$edge)

  mean_mu <- list()
  mean_lambda <- list()
  for (i in 1:n_edges){
    Em <- dplyr::select(x$forward[[i]], starts_with("E")) %>% as.matrix()
    Dm <- dplyr::select(x$forward[[i]], starts_with("D")) %>% as.matrix()
    Fm <- dplyr::select(x$forward[[i]], starts_with("F")) %>% as.matrix()

    Pm = (Fm * Dm) / (rowSums(Fm * Dm))

    mean_mu[[i]] <- sum(rowSums(sapply(1:k, function(i) Pm[,i] * mu[i]))) / ntimes
    mean_lambda[[i]] <- sum(rowSums(sapply(1:k, function(i) Pm[,i] * lambda[i]))) / ntimes
  }

  mean_rates <- list(
    "mu" = unlist(mean_mu),
    "lambda" = unlist(mean_lambda)
  )

  res <- list(
    "mean_rates" = mean_rates,
    "x" = x
  )

  return(res)
}
