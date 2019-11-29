\name{gen.factorial}
\alias{gen.factorial}
\title{Generates a full factorial design}
\description{
 A full factorial design is generated.
}
\usage{
	gen.factorial(levels, nVars=0, center=TRUE, factors="none",varNames=NULL)
}

\arguments{
	\item{levels}{A vector of levels for the variables. May be an integer if nVars is specified.}
	\item{nVars}{The number of variables.}
	\item{center}{If TRUE, all non-factors will be centered.}
	\item{factors}{If "all", all variables are factors, otherwise a vector of the variable numbers
		of the variables that are to be factors.}
	\item{varNames}{The names of the variables.}
}
\value{
	A factorial design in a data.frame.
}

\examples{
dat<-gen.factorial(3,3)
dat<-gen.factorial(c(3,2,3))
dat<-gen.factorial(3,3,factors="all")
dat<-gen.factorial(3,3,varNames=c("A","B","C"))

}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). gen.factorial. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}


\keyword{design}

