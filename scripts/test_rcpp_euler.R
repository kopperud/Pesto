library(microbenchmark)
library(BDS)


microbenchmark(
 BDS:::rcpp_backwards(lambda, mu, eta, c(0.0, 0.0, 1.0, 1.0), 21.1235, 1000)
)

l <- list()
step_sizes <- floor(slouch::lseq(from = 10, to = 100000, length.out = 15))

for (i in seq_along(step_sizes)){
  nsteps <- step_sizes[i]

  u0 <- c(0.0, 0.0, 1.0, 1.0)
  bl <- 21.5
  l[[i]] <- BDS:::rcpp_backwards(lambda, mu, eta, bl, nsteps)
}

