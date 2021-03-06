\name{RCLM}
\alias{RCLM}
\title{ Fit Random Coefficient Linear Models}
\description{An extension of the ICRO algorithm for Bayesian Computation. It can be used to fit a Random Coefficient Linear Models and estimate the coefficients \eqn{\beta} and \eqn{\sigma^2}.}

\usage{
RCLM(I=100, J=10, Data, iteration = 10000, warm = 100)
}

\arguments{
  \item{ I }{ Number of first subjects in the random coefficient linear model (RCLM).}
  \item{ J }{ Number of second subjects in the random coefficient linear model (RCLM).}
  \item{ Data }{ A simulated dataset. The first column is the response and the rest is for explanatory variables, see \code{RCDat} for detail.}
  \item{ iteration  }{ The number of total iterations, the default value is 10000.}
\item{ warm  }{ The number of burn-in iterations, the default value is 100.}

}

\value{
  \item{path}{The traces of estimated coefficients vs. iterations.}
  \item{coef}{The mean of estimated coefficients \eqn{\mathbf{\beta}} and \eqn{\sigma^2}.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
\donttest{
library(IROmiss)
data(RCDat)
RCLM(I=100, J=10, RCDat, iteration = 10000, warm = 1000)
}

         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Penalized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{RCLM}
