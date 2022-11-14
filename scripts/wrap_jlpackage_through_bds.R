## Install JuliaCall, Julia and the other package
install.packages("JuliaCall")
library(JuliaCall)

install_julia()
julia_setup()
julia_eval('import Pkg; Pkg.add(url = "https://github.com/kopperud/Diversification.jl.git")')

## Install the R-package
## (will need a working SSH-key since the repository is private)
remotes::install_git("git@github.com:hoehna/BDS_deterministic_map.git", subdir = "BDS")

## Run the analysis using the R-package wrapper
library(BDS)
library(ape)

data(primates)

lambda <- c(0.1, 0.2)
mu <- c(0.05, 0.15)
rho <- 1.0
eta <- 0.05

res <- birth_death_shift3(primates, lambda, mu, eta, rho, verbose = TRUE)
lambda_average <- res$lambda
