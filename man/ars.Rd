\name{ars}
\alias{ars}
\title{Adaptive Rejection Sampling for Gibbs Sampling}
\usage{
ars(g, n, lower = -Inf, upper = Inf)
}
\description{
ars is used to conduct adaptive rejection sampling on a valid probability density g which is log concave. 
}
\arguments{
\item{g}{A function denoting a valid probability density on some interval defined by the lower and upper arguments. The function will check if it is log concave on the interval, and if not will produce an error message}
\item{n}{The number of observations to be drawn from probability density g.}
\item{lower}{The lower bound of the interval to sample from density g on. The default is -Inf, and if the user doesn't specify bounds will be autmatically determined by the function.}
\item{upper}{The upper bound of the interval to sample from density g on. The default is Inf, and if the user doesn't specify bounds will be autmatically determined by the function.}
}

\examples{
X1 <- ars(g = dnorm, n = 1e3, lower = -5, upper = 5)      
X2 <- ars(g = dexp, n = 1e3, lower = 0, upper = 10)

## will produce an error since Cauchy distribution is not log concave
X3 <- ars(g = dcauchy, n = 1e3, lower = 0, upper = 20)

## sets default bounds to -38 and 38 since unspecified
X4 <- ars(g = dnorm, n = 1e3)
}
