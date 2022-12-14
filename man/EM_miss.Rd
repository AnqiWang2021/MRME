% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM_miss.R
\name{EM_miss}
\alias{EM_miss}
\title{Estimate parameters using the expectation-maximization algorithm}
\usage{
EM_miss(para_start, Sigma_start, G, Y, X, maxit, tol)
}
\arguments{
\item{para_start}{numeric vectors of starting values for beta and gamma. Beta is the parameter in the linear regression model of
response variable and exposure variable of interest. Gamma is the parameter in the linear regression model of exposure variable
and genetic matrix.}

\item{Sigma_start}{numeric vectors of starting values for Sigma, Sigma is the covariance matrix of the error terms.}

\item{G}{A numeric matrix of genetic variants}

\item{Y}{A numeric vector of the outcome variable}

\item{X}{A numeric vector of the exposure variable}

\item{maxit}{Maximum number of iterations to solve the estimating equations}

\item{tol}{Numerical precision}
}
\value{
A list
\describe{
\item{step}{number of iterations of the estimation algorithm}
\item{Sigma}{the final estimates of the covariance matrix.}
\item{para}{the final estimates of the parameters organized as (beta and gamma).}
}
}
\description{
\code{EM_miss} Estimate parameters using the expectation-maximization algorithm based on the complete log-likelihood
of observed exposure data in the presence of missing data.
}
\examples{
library(MRME)
library(MASS)
N = 10000  #the sample size
K = 1     #the number of instrumental variables
n = 2
p = 0.3
#Sigma matrix
Sigma = matrix(c(0.58,0.35,0.35,0.58),2,2)
beta0 = 0.8
beta1 = 0.2
gamma0 = 0.13
gammav = 0.2
epsilon = mvrnorm(N,rep(0,2),Sigma)
G<-rbinom(N*K,n,p)
G = as.matrix(G)
X = gamma0 + as.matrix(G)\%*\%gammav + epsilon[,1]
Y = beta0 + beta1 *X +epsilon[,2]
p_miss =  exp(0.25 + 0.6*Y +  G\%*\%rep(0.05,1))/(1 + exp(0.25 +0.6*Y + G\%*\%rep(0.05,1)))
maxit = 1000
tol = 1e-8
Sigma_start= matrix(c(0.4,0.2,0.2,0.4),2,2)
para_start = c(0.6,0.1,rep(0.05,2))
result = EM_miss(para_start,Sigma_start,G,Y,X,maxit,tol)
}
