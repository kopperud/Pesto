% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_eta.R
\name{optimize_eta}
\alias{optimize_eta}
\title{Estimate the shift rate conditional on lambda and mu}
\usage{
optimize_eta(lambda, mu, phy, rho, lower = 1e-05, upper = 10, eta0 = NULL)
}
\arguments{
\item{lambda}{a vector of speciation rates}

\item{mu}{a vector of extinction rates}

\item{phy}{an object of class phylo}

\item{lower}{lower limit for the numerical optimization of eta}

\item{upper}{upper limit for the numerical optimization of eta}

\item{eta0}{starting value for eta}
}
\value{

}
\description{
finds the maximum-likelihood of eta conditional on the vectors lambda and mu, using gradient descent. This will typically evaluate the likelihood function around 20-60 times, and also the derivative of the likelihood function with respect to eta, so it is not instantaneous. Check if the result is stuck at the lower or upper bounds.
}
\examples{
library(Pesto)
data(primates)

lambda <- c(0.1, 0.2, 0.3, 0.4)
mu <- c(0.05, 0.15, 0.05, 0.15)
rho <- 0.67

etaml <- optimize_eta(lambda, mu, primates, rho)
print(etaml)
}
