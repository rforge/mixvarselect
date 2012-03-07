\name{moeridge}
\alias{moeridge}
\title{
mixture of expert fit with ridge penalization
}
\description{
This function fits a two-component mixture of expert model with ridge penalization on parameters of the linear predictor of proportions.
}
\usage{
moeridge(y, x, lambda.ridge = 0.1 * log(nrow(x)), 
initial = list(alpha = rep(0.5, (ncol(x) + 1)), 
beta1 = rep(0, ncol(x)), beta2 = rep(0, ncol(x)), 
sigma = 0.88), conv.eps = 1e-08, 
maxiter = list(total = 2500, em = 15))
}
\arguments{
  \item{y}{
numeric vector of response to fit.}
  \item{x}{
matrix of predictors with samples in rows and covariates in columns.
}
  \item{lambda.ridge}{
ridge penalization constant.}
  \item{initial}{
a list of initial values.}
  \item{conv.eps}{
EM convergence threshold.}
  \item{maxiter}{
a list of prespecified integer values for maximum number of iterations.}
}
\details{
we may fill this later
}
\value{
\item{alpha}{the parameters of the linear predictor of proportions under logit transformation.}
\item{beta}{matrix of parameters of the linear predictors for mean components.}
}
\references{
Khalili, A. (2010). New Estimation and Feature Selection Methods in Mixture-of-Experts Models. The Canadian Journal of Statistics, 38, 519-539.
}
\author{
Vahid Partovi Nia and Abbas Khalili
}
\note{
The choice of initial values is crucial. A common variance is assumed. Maximized log likelihood value may be returned, or the estimated common variance. This function only fits Gaussians with two components. 
A good choice of penalization constant is proved to be of order log of sample size. This function requires installation of GSL libraries on the machine and may not work under Windows platform.
}

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

\examples{
data(moesim)
moeridge(moesim$y,moesim$x)
}
\keyword{mixture of expert}
\keyword{ridge penalty}
