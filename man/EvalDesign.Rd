\name{eval.design}
\alias{eval.design}
\title{Evaluates a design.}
\description{
 A design is evaluated.
}
\usage{
eval.design(frml,design,confounding=FALSE,variances=TRUE,center=FALSE,X=NULL)
}
\arguments{
	\item{frml}{The formula used to create the design.}
	\item{design}{The design, which may be the design part of the output of optFederov().}
	\item{confounding}{If confounding=TRUE, the confounding patterns will be shown.} 
	\item{variances}{If TRUE, the variances each term will be output.}
	\item{center}{If TRUE, numeric variables will be centered before frml is applied.}
	\item{X}{X is either the matrix describing the prediction space for I or for G, the 
		the candidate set from which the design was chosen. They are often
		the same.}
}
\note{
	I, Ge and Dea are calculated only when X is input.
}
\value{
	\item{confounding}{A matrix. The columns of which give the regression coefficients of
		each variable regressed on the others. If \eqn{C} is the confounding matrix, then 
		\eqn{-ZC} is a matrix of residuals of the variables regressed on the other variables.}
	\item{determinant}{\eqn{|M|^{1/k}}{(det(M/N)^(1/k)}, where \eqn{M=Z'Z/N}{M=Z'Z/N}, and Z is 
		the model expanded \eqn{N\times k}{N x k} design matrix.}
  	\item{A}{The average coefficient variance: \eqn{trace(Mi)/k}, where \eqn{Mi} is the 
		inverse of \eqn{M}.}
  	\item{I}{The average prediction variance over X, which can be shown to be 
		\eqn{trace((X'X*Mi)/N)}{trace((X'X*Mi)/N}.} 
  	\item{Ge}{The minimax normalized variance over X, expressed as an efficiency with respect 
		to the optimal approximate theory design. It is defined as \eqn{k/max(d)}, where 
		\eqn{max(d)} is the maximum normalized variance over \eqn{X} -- i.e. the max of 
		\eqn{x'(Mi)x}, over all rows \eqn{x'} of \eqn{X}.}
  	\item{Dea}{A lower bound on \code{D} efficiency for approximate theory designs. It is 
		equal to \eqn{exp(1-1/Ge)}.}
	\item{diagonality}{The diagonality of the design, excluding the constant, if any. Diagonality
		is defined as \eqn{(|M_1|/\prod{diag(M_1)})^{1/k}}{(|M1|/prod(diag(M1)))^(1/k)}, where 
		\eqn{M_1}{M1} is \eqn{M} with first column and row deleted when there is a constant.}
	\item{gmean.variances}{The geometric mean of the coefficient variances.}

}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). eval.design. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}

\keyword{design}
