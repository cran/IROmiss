\name{RCDat}
\alias{RCDat}
\docType{data}

\title{A simulated dataset for random coefficient models.} 

\description{
The dataset is generated using the default settings. The Number of customers I=100 and each customer responds to J=10 items. For the parameters, the true coefficient \eqn{\mathbf{\beta}} is \eqn{(\beta_0,\beta_1,\beta_2,\beta_3)= (1, 2, 1.5, 1)} and the true value of \eqn{\sigma^2} is 0.25. The first column of the dataset denote the response \eqn{\mathbf{y}}. The dataset should be used in \code{RCLM(RCDat)}.
\describe{
  \item{ RCDat }{ A simulated dataset.}
}
}
\usage{data(RCDat)}
\format{matrix}




\references{
Liang, F., Jia, B., Xue, J., Li, Q., and Luo, Y. (2017). An Imputation Regularized Optimization Algorithm for High-Dimensional Missing Data Problems and Beyond. Submitted to Journal of the Royal Statistical Society Series B.
}


\keyword{datasets}
