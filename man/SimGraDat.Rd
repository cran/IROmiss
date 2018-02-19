\name{SimGraDat}
\alias{SimGraDat}
\title{ Simulate Incomplete Data for Gaussian Graphical Models}
\description{Simulate incomplete data with a band structure, which can be used in \code{GraphIRO(data,...)} for estimating the structure of the Gaussian graphical network.}

\usage{
SimGraDat(n = 200, p = 100, type = "band", rate = 0.1)
}

\arguments{
  \item{ n }{ Number of observations, default of 200.}
  \item{ p }{ Number of covariates, default of 100.}
  \item{ type }{\code{type=="band"} which denotes the band structure, with precision matrix
  \deqn{
        C_{i,j}=\left\{\begin{array}{ll}
                     0.5,&\textrm{if $\left| j-i \right|=1, i=2,...,(p-1),$}\\
                     0.25,&\textrm{if $\left| j-i \right|=2, i=3,...,(p-2),$}\\
                     1,&\textrm{if $i=j, i=1,...,p,$}\\
                    0,&\textrm{otherwise.}
                \end{array}\right. }
                }
  \item{ rate }{ Missing rate, the default value is 0.1.}
}

\value{
   \item{ data }{ \emph{n}x\emph{p} Gaussian distributed data with missing.}
  \item{ A }{  \emph{p}x\emph{p} adjacency matrix used for generating data.}
%% ...
}

\author{ Bochao Jia\email{jbc409@ufl.edu} and Faming Liang}

\examples{
library(IROmiss)
SimGraDat(n = 200, p = 100, type = "band", rate = 0.1)
}

\references{
Liang, F., Song, Q. and Qiu, P. (2015). An Equivalent Measure of Partial Correlation Coefficients for High Dimensional Gaussian Graphical Models. J. Amer. Statist. Assoc., 110, 1248-1265.

Liang, F. and Zhang, J. (2008) Estimating FDR under general dependence using stochastic approximation.  Biometrika,  95(4), 961-977.

Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2018). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.

}

\keyword{SimGraDat}
