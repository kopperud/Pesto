library(microbenchmark)
library(BDS)
library(tibble)
library(ggplot2)

## How many nsteps to get good logL?

data(primates)
data(bears)

lambda <- c(0.1, 0.2)
mu <- c(0.05, 0.15)
eta <- 0.05

datasets <- list("primates" = primates,
             "bears" = bears)

jl_logL <- list("bears" = -20.628120241987435,
                "primates" = -696.2978988828179)
dfs <- list()
for (dataset_name in names(datasets)){
  phy <- datasets[[dataset_name]]

  step_sizes <- floor(seq(100, 5000, length.out = 20))
  logLs <- list()
  times <- list()
  for (i in seq_along(step_sizes)){
    cat(".",i)
    nsteps <- step_sizes[i]
    time1 <- Sys.time()
    res <- birth_death_shift2(phy, lambda, mu, eta, ntimes = nsteps)
    time2 <- Sys.time()
    times[[i]] <- time2 - time1
    logL <- res[["logL"]]
    logLs[[i]] <- logL
  }; times <- unlist(times)


  df1 <- tibble("logL" = unlist(logLs),
                "nsteps" = step_sizes,
                "times" = times,
                "implementation" = "Rcpp",
                "dataset" = dataset_name)
  df2 <- tibble("logL" = rep(jl_logL[[dataset_name]], 2),
                "nsteps" = range(step_sizes),
                "times" = NaN,
                "implementation" = "Julia",
                "dataset" = dataset_name)
  df <- dplyr::bind_rows(df1, df2)
  dfs[[dataset_name]] <- df
}
df <- dplyr::bind_rows(dfs)

p1 <- ggplot(df, aes(x = nsteps, y = logL, col = implementation)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = c(0.2, 0.5)) +
  xlab("Number of steps in Euler's method \n(per branch, dt = branch_length / n)") +
  scale_color_manual(values = c("orange", "black")) +
  facet_wrap(~ dataset, scales = "free")


p2 <- ggplot(df, aes(x = nsteps, y = times)) +
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = c(0.2, 0.5)) +
  xlab("Number of steps in Euler's method \n(per branch, dt = branch_length / n)") +
  scale_color_manual(values = c("black", "orange")) +
  facet_wrap(~ dataset, scales = "free") +
  ylab("runtime (seconds)")

library(patchwork)

p <- p1 / p2

ggsave("../figures/rcpp_euler_vs_jl.pdf", units = "mm", height = 200, width = 200)


plot(step_sizes, unlist(logLs))
abline(h =  -696.2978988828179, col = "red")






l <- list()
step_sizes <- floor(slouch::lseq(from = 10, to = 100000, length.out = 15))
times <- list()

for (i in seq_along(step_sizes)){
  nsteps <- step_sizes[i]

  time1 <- Sys.time()
  l[[i]] <- BDS:::birth_death_shift2(primates, lambda, mu, eta, ntimes = nsteps)
  time1 <- Sys.time()
  times[[i]] <- time2 - time1
}

##




