\name{model.matrix.formula}
\alias{model.matrix.formula}
\title{Builds a model matrix}
\description{
 Produces a model matrix using expand.formula()
}
\usage{
\method{model.matrix}{formula}(frml,data,...)
}

\arguments{
	\item{frml}{A formula.}
	\item{data}{A data.frame}
	\item{...}{Additional arguments passed to model.matrix.default()}
}
\details{
	This is a method function signaled by a formula as the first argument. 
	It causes the formula to be translated by expand.formula() before calling
      model.matrix.default().
}
\value{
	A matix.
}
\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). model.matrix.formula. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}

\keyword{design}

