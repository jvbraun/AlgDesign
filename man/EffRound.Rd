\name{efficient.rounding}
\alias{efficient.rounding}
\title{Efficient Rounding}
\description{
 A vector of proportions is efficiently rounded to integers.
}
\usage{
efficient.rounding(proportions,n,random=TRUE)
}

\arguments{
	\item{proportions}{A vector of proportions.}
	\item{n}{The sum of the resulting integers.}
	\item{random}{If TRUE, ties will be broken at random, otherwise the first of the tied values will be used.}
}
\details{
	This function implements an efficient rounding procedure to round approximate theory designs
	into replicated integer approximations.
}
\value{
	A vector of replications summing to n.
}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). efficient.rounding. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}


\references{
  Pulkesheim, F. and Rieder, S. (1992). Efficient rounding of approximate designs. Biometrika. 79-4. 763-770.}
\keyword{design}
