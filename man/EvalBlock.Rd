\name{eval.blockdesign}
\alias{eval.blockdesign}
\title{Evaluates a blocked design.}
\description{
 A blocked design is evaluated.
}
\usage{
eval.blockdesign(frml,design,blocksizes,rho=1,confounding=FALSE,center=FALSE)
}
\arguments{
	\item{frml}{The formula used to create the blocked design.}
	\item{design}{The blocked design, which may be the design output by optBlock().}
	\item{blocksizes}{A vector of blocksizes for the design.}
	\item{rho}{A vector, giving the ratios of whole to within variance components.}
	\item{confounding}{If confounding=TRUE, the confounding matrix will be output.} 
	\item{center}{If TRUE, numeric variables will be centered before frml is applied.}

}
\value{
	\item{confounding}{A matrix based on the design matrix in which the within block variables
		have been centered about their block means. The columns of this matrix which give 
		the regression coefficients of each variable regressed on the others. If \eqn{C} is the 
		confounding matrix, then \eqn{-XC} is a matrix of residuals of the variables regressed on 
		the other variables.}
	\item{determinant.all.terms.within.terms.centered}{\eqn{|M|^{1/k}}{(det(M)^(1/k)}, where 
		\eqn{M=X^TX/N}{M=X'X/N} and X is the model expanded \eqn{N\times k}{N x k} design matrix 
		in which the within block variables	have been centered about the grand mean.}
	\item{within.block.efficiencies}{The determinant criterion blocking efficiencies for the range 
		of rho's input. A high efficiency indicates that there is little intrablock information
		to be recovered in the analysis.}
	\item{block.centered.properties}{A matrix with four rows. The columns correspond to constant, whole
		block terms and within block terms:

		\enumerate{
			\item The degrees of freedom for terms in the expanded model.  
			\item The determinant of the block centered within block terms.  
			\item The geometric mean of the block centered variances.  
			\item The geometric mean of the ratio of centered to block centered variances. 
			}
		}
}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). eval.blockdesign. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}

\keyword{design}

