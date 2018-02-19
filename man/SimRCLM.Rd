\name{SimRCLM}
\alias{SimRCLM}
\title{ Simulate Dataset for Random Coefficient Linear Models}
\description{Simulate a dataset for random coefficient linear model, which can be used in \code{RCLM(I,J,RCDat,...)}.}

\usage{
SimRCLM(I=100, J=10, beta, sigma)
}

\arguments{
  \item{ I }{ Number of first subjects in the random coefficient linear model (RCLM).}
  \item{ J }{ Number of second subjects in the random coefficient linear model (RCLM).}
  \item{ beta }{ A 4x1 vector of random coefficients of the model, now only allows length 4.}
  \item{ sigma  }{ The standard diviation for the noise term.}

}

\value{
  \item{D}{A simulated data matrix for random coefficient models. The first column of the dataset denote the response \eqn{\mathbf{y}}. The dataset should be used in \code{RCLM(I,J,RCDat...)}.}
  \item{coef}{The mean of estimated coefficients \eqn{\mathbf{\beta}} and \eqn{\sigma^2}.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
\dontshow{
library(IROmiss)
beta<-c(1,-1,3,-1.5)
sigma <- 0.5
D <- SimRCLM(I=75, J=15, beta, sigma)
D
RCLM(I=75, J=15, D, iteration = 100, warm = 10)
}
\donttest{
library(IROmiss)
beta<-c(1,2,1.5,1)
sigma <- 0.5
D <- SimRCLM(I=100, J=10, beta, sigma)
RCLM(I=100, J=10, D, iteration = 10000, warm = 1000)
}

}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Penalized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{RCLM}
