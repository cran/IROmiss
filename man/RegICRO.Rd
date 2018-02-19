\name{RegICRO}
\alias{RegICRO}
\title{ Variable selection for high-dimensional Regression with Missing Data.}
\description{Application of the imputation conditional regularized optimization (ICRO) algorithm for high-dimensional variable selection in presence of missing data.}

\usage{
RegICRO(x, y, coef, type = "indep", alpha1 = 0.1, alpha2 = 0.05,
iteration = 30, warm = 20)
}

\arguments{
  \item{ x }{  \emph{\eqn{n}}x\emph{p} covariates matrix.}
  \item{ y }{  \emph{\eqn{n}}x1 responses.}
  \item{ coef }{  A \emph{p}x1 vector of coefficients for the linear regression model. The intercept coefficient is default to 1.}
    \item{ type }{When \code{type=="indep"}, the case with independent covariates, or \code{type=="dep"}, the case with dependent covariates, the default type is "indep".}
  \item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see R package \pkg{equSA} for detail. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.1.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.05.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm  }{ The number of burn-in iterations, the default value is 20.}
}

\value{
\item{Var}{Selected variables and their estimated coefficients by our ICRO algorithm.}
  \item{table}{The summarized table for evaluating the performance of ICRO algorithm. 'bias' denotes Euclidean distance between estimated coefficients and true coefficients; 'fsr' denotes false selection rate and 'nsr' denotes negative selection rate. The smaller the measurements are, the better the performance is.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
      \dontshow{
library(IROmiss)
p <- 50
beta <- rep(0,p)
beta[1:5] <- c(1, 2, -1.5, -2.5, 5)
result <- SimRegDat(n = 100, p = 50, coef = beta, data.type = "indep",
miss.type="MAR", rate = 0.05)
RegICRO(result$x, result$y, result$coef, type = "indep", iteration = 10, warm = 5)
      }
    \donttest{
library(IROmiss)
p <- 200
beta <- rep(0,p)
beta[1:5] <- c(1, 2, -1.5, -2.5, 5)
result <- SimRegDat(n = 100, p = 200, coef = beta, data.type = "indep", 
miss.type="MAR", rate = 0.05)
RegICRO(result$x, result$y, result$coef, type = "indep", iteration = 20, warm = 10)
      }
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{RegICRO}
