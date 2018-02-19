\name{eye_norm}
\alias{eye_norm}
\docType{data}

\title{Example dataset for high-dimensional variable selection by the ICRO algorithm.} 

\description{
Normalized Gene expression data from the microarray experiments of mammalian-eye tissue samples of Scheetz et al. (2006). It should be used in \code{EyeICRO(x,y...)}.
\describe{
  \item{ x }{ a \emph{n}x\emph{p} gene expression data.}
  \item{ y }{The expression level of gene TRIM32.}
}
}
\usage{data(eye_norm)}
\format{A list containing the matrix \code{x} and response matrix
  \code{y}}




\references{
T. Scheetz, k. Kim, R. Swiderski, A. Philp, T. Braun, K. Knudtson, A. Dorrance, G. DiBona, J. Huang, T. Casavant, V. Sheffield, E. Stone .Regulation of gene expression in the mammalian eye and its relevance to eye disease. Proceedings of the National Academy of Sciences of the United States of America, 2006.
}


\keyword{datasets}
