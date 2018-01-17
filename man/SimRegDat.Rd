\name{SimRegDat}
\alias{SimRegDat}
\title{ Simulate Incomplete Data for High-Dimensional Linear Regression.}
\description{Simulate incomplete data for high-dimensional linear regression with dependent or independent covariates\code{RegICRO(x,y...)}.}

\usage{
SimRegDat(n = 100, p = 200, coef, data.type = "indep",
miss.type="MCAR", rate = 0.1)
}

\arguments{
  \item{ n }{ Number of observations, default of 100.}
  \item{ p }{ Number of covariates, default of 200.}
  \item{ coef }{ A \emph{p}x1 vector of coefficients for the linear regression model. The intercept coefficient is default to 1.}
  \item{ data.type }{ When \code{data.type=="indep"}, it simulates the data with independent covariates, each covariate independently follow the normal distribution with mean 0 and variance 4. When \code{data.type=="dep"}, it simulates the data with dependent covariates with "band" dependent structure, see \code{SimGraDat} for detail. The default data type is "indep".}
  \item{ miss.type }{ \code{miss.type=="MCAR"} refer to the case of missing completely at random. when \code{miss.type=="MAR"}, the missing probability for each data point is proportional to the mean of its conditional normal distribution, the default missing type is "MCAR".}
  \item{ rate }{ Missing rate, the default value is 0.1.}
}

\value{
   \item{ x }{ \emph{n}x\emph{p} covariates matrix.}
  \item{ y }{  \emph{n}x1 responses.}
  \item{ coef }{ \emph{p}x1 vector of coefficients for the linear regression model.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
library(IROmiss)
p <- 200
beta <- rep(0,p)
beta[1:5] <- c(1, 2, -1.5, -2.5, 5)
SimRegDat(n = 100, p = 200, coef = beta, data.type = "dep", 
miss.type="MAR", rate = 0.1)
         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{SimRegDat}
