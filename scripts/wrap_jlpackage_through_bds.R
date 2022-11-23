## Install JuliaCall, Julia and the other package
install.packages("JuliaCall")
library(JuliaCall)

install_julia()
julia_setup()
julia_eval('import Pkg; Pkg.add(url = "https://github.com/kopperud/Diversification.jl.git")')

## Install the R-package
## (will need a working SSH-key since the repository is private)
## remotes::install_git("git@github.com:hoehna/BDS_deterministic_map.git", subdir = "BDS")
remotes::install_local("~/projects/BDS_deterministic_map/BDS")

## Run the analysis using the R-package wrapper
library(BDS)
library(ape)

data(primates)

lambda <- c(0.1, 0.2)
mu <- c(0.05, 0.15)
rho <- 1.0
eta <- 0.05

microbenchmark::microbenchmark(birth_death_shift3(primates, lambda, mu, eta, rho, verbose = F))
lambda_average <- res$lambda


library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)


th <- max(node.depth.edgelength(primates))

df1 <- tibble("node" = 1:max(primates$edge),
              "Speciation rate" = lambda_average)
x <- as_tibble(primates)

phydf <- merge(x, df1, by = "node")
td_phy <- as.treedata(phydf)

p1a <- ggtree(td_phy, aes(color = `Speciation rate`)) +
  geom_tiplab(size = 8) +
  theme(legend.position = c(0.2, 0.8)) +
  xlim(c(0.0, th + 10))

##
ggsave("/tmp/p1a.pdf", p1a)
