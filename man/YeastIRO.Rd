\name{YeastIRO}
\alias{YeastIRO}
\title{ Learning gene regulatory networks for Yeast Cell Expression Data.}
\description{An Imputation Regularized Optimization (IRO) algorithm for learning gene regulatory networks with missing data. The dataset is collected from the yeast Saccharomyces cerevisiae responding to diverse environmental changes and is available at http://genome-www.stanford.edu/yeast-stress/.}

\usage{
YeastIRO(data, alpha1 = 0.05, alpha2 = 0.01, alpha3 = 0.01, iteration = 30, warm = 20)
}

\arguments{
  \item{ data }{ \emph{\eqn{n}}x\emph{p} Yeast Cell expression data.}
 \item{ alpha1  }{ The significance level of correlation screening in the \eqn{\psi}-learning algorithm, see R package \pkg{equSA} for detail. In general, a high significance level of correlation screening will lead to 
 a slightly large separator set, which reduces the risk of missing important variables in 
 the conditioning set. In general, including a few false variables in the conditioning set will not hurt much the 
 accuracy of the \eqn{\psi}-partial correlation coefficient, the default value is 0.05.}
  \item{ alpha2  }{ The significance level of \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix, see \pkg{equSA}, the default value is 0.01.}
  \item{ alpha3  }{ The significance level of integrative \eqn{\psi}-partial correlation coefficient screening for estimating the adjacency matrix of IRO_Ave method, the default value is 0.01.}
\item{ iteration  }{ The number of total iterations, the default value is 30.}
\item{ warm  }{ The number of burn-in iterations, the default value is 20.}
}

\value{
  \item{ A }{  \emph{p}x\emph{p} Estimated adjacency matrix for network construction.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
\dontshow{
library(IROmiss)
library(huge)
data(yeast)
huge.plot(diag(100))
}
\donttest{
library(IROmiss)
library(huge)
data(yeast)
## long time ##
A <- YeastIRO(yeast, alpha1 = 0.05, alpha2 = 0.01, alpha3 = 0.01, iteration = 20, warm = 10)
## plot gene regulatory network by our estimated adjacency matrix.
huge.plot(A)
}

         
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{YeastIRO}
