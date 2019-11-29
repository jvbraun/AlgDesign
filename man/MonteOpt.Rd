\name{optMonteCarlo}
\alias{optMonteCarlo}
\title{Optimal design via Monte Carlo}
\description{
  Finds a design using the specified criterion via Federov's algorithm
  applied to a random subset of all possible candidate points. 
}
\usage{
optMonteCarlo(frml,data,nTrials,approximate=FALSE,criterion="D",evaluateI=FALSE,
	space=NULL,mixtureSum=1,constraints=NULL,RandomStart=TRUE,nRepeats=5,nCand,
	nCandNull,DFrac=1,CFrac=1,args=FALSE)
}


\arguments{
	\item{frml}{Required: A formula starting with ~ which will be used with 
		model.matrix() to create a model matrix. If there are mixture variables,
		the constant term is suppressed.}
	\item{data}{Required: A data frame with 7 or 8 columns. See details below
		for specifics} 
	\item{nTrials}{number trials in design -- must be greater than the
	  number of terms in the model, if  missing will be set to the
	  number of model terms in the model  plus five.}
	\item{approximate}{When FALSE, an exact design in nTrails will be calculated. When TRUE the proportions for
	 an approximate theory design will be calculated. If nTrials is set, any proportion less than 1/(2*maxIteration)
	 will be discarded before the proportions are efficiently rounded, otherwise all non-zero proportions will be 
	 shown: these are the support points.}
	\item{criterion}{"D", "A", or "I"}
	\item{evaluateI}{TRUE if I is to be evaluated in addition to the
	  other criteria -- slower because of calculations for I}
	\item{space}{If the criterion is "I" or evaluate I is true, the space over which the I criterion is to
	  be evaluated may be input. It should be a matrix with the same column types and names as in data. If
	  space is not input the evaluation will be done over the space described by data.}	
	\item{constraints}{A function taking a vector argument with length
	  equal to the number of variables, and returning TRUE if the vector is 
	  inside the constrained region}
	\item{mixtureSum}{The mixture variables, if any, will sum to this value.}
	\item{RandomStart}{When TRUE, the starting design will be chosen
	  at random, otherwise nullification will be used. Note: the nullifcation
	  used here is different and much slower than that in optFederov().}
	\item{nRepeats}{number of times to retry the entire process}
	\item{nCand}{number of candidate points to generate, if missing,
	  it will be 10 times the number of terms}
	\item{nCandNull}{Number of candidate points to use for
	  nullification. If missing it will be set to nCand}
	\item{DFrac}{Fraction of design used in search: 1 uses all of
	  them, 0 only the one with the smallest variance}
	\item{CFrac}{Fraction of candidate set searched : 1 uses all of
	  them, 0 only the one with the largest variance}
	\item{args}{If TRUE, the actual arguments to the function including the starting 
		random number seed will be output.}
}
\details{
  The columns of the input data frame are as follows. The columns need
	not be named. It is probably best to avoid naming the variables with
	single letters, especially "I" -- use paste(), as in the examples. For each variable
	nLevels are randomly generated between low and high, inclusive, and then rounded
	with round. For integer levels, round should be set to 0.

	\describe{
	\item{var:}{The names of the variables.}
	\item{low:}{The lower limit of the range for each variable. Ignored for mixtures.} 
	\item{high:}{The upper limit of the range for each variable. Ignored for mixtures.} 
	\item{center:}{The centering value for each variable. Ignored for mixtures.}
	\item{nLevels:}{The number of levels for each variable. Ignored for mixture variables.}
	\item{round:}{The number of decimal digits for the levels. The levels are randomly and 
		uniformly chosen between low and high, and this parameter controls the number of 
		trailing digits. The max value for mixture variables in this vector is used to 
		round all mixture variables.}
	\item{factor:}{TRUE, FALSE depending on whether or not the variable is a factor. Note: other
		columns will be reset to conform to a nLevels factor.}
	\item{mix:}{TRUE if the variable is a mixture variable. This column may
		be omitted if there are no mixture variables.}
	}

  Candidate lists required by \code{optFederov()} increase with the number of
  variables, and can easily exceed storage capacity and can require
  excessive amounts of time to process. To overcome this problem,
  \code{optMonteCarlo()}, generates at random  \code{nCand}  points from a putative 
  candidate list.

  For non-mixture variables, \code{optMonteCarlo()} samples from the putative
  candidate list by choosing random levels inside the limits given by
  \code{low} and \code{high} in \code{data}. These are rounded to the 
  number of levels given by \code{nLevels} in \code{data} and to the 
  number of decimal digits given by \code{round} in \code{data}.

  For mixture variables, \code{optMonteCarlo()} samples from the putative
  candidate list by choosing random levels between 0 and 1, rounded to
  the maximum in the \code{round} column of \code{data}, and such that the sum over 
  all variables is equal to \code{mixtureSum}.
  
  If a constraint function is supplied in \code{Constraints}, it is
  applied, and results which do not meet the constraint are
  discarded. The constraint function should be written to process 
  uncentered variables.

  The above procedures are repeated until \code{nCand} candidate points
  are found.

  Nullification, successively adds points to a design until n points are 
  found. This is the same procedure that is in \code{optFederov} except
  that each new point is selected from a new sampling of the putative
  candidate points. In general, this will produce better designs that
  those from a random start.

  The entire process is repeated \code{nRepeats} times, and the best
  result is reported. The methodology compares favorably with an
  exhaustive search where the entire candidate list is searched by
  \code{optFederov()}.
  
  The random numbers used in these calculations are controlled by the
  usual R random number mechanism.

	A vignette giving further details is availble. To access it, type
	
	vignette("AlgDesign")

}



\value{
The output is the same list as from \code{optFederov}, but the criteria
values are relative to the randomly chosen subsets of the putative candidate 
space. In general, they should not differ greatly from those obtained by an 
exhaustive search.

  \item{D}{The kth root of the generalized variance: \eqn{det(M)^{1/k}}{det(M)^(1/k)}, where \eqn{det(M)} is the
    determinant of the normalized dispersion matrix \eqn{M} --
    i.e. \eqn{M=Z'Z/n}, where \eqn{Z=X[rows,]}}
  \item{A}{The average coefficient variance: \eqn{trace(Mi)/k}, where \eqn{Mi} is the inverse of \eqn{M}.}
  \item{I}{The average prediction variance over X, which can be shown to be \eqn{trace((X'X*Mi)/N)}. 
	This is calculated only when I is the criterion or when \code{evaluateI} is TRUE.}
  \item{Ge}{The minimax normalized variance over X, expressed as an efficiency with respect to the optimal approximate
	theory design. It is defined as \eqn{k/max(d)}, where \eqn{max(d)} is the maximum normalized
    	variance over \eqn{X} -- i.e. the max of \eqn{x'(Mi)x}, over all rows \eqn{x'}
    	of \eqn{X}.}
  \item{Dea}{A lower bound on \code{D} efficiency for approximate theory designs. It is equal to \eqn{exp(1-1/Ge)}.}
  \item{Design}{The design.} 
  \item{args}{A list of the actual arguments used in this call.}

}
\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). optMonteCarlo. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}

\examples{

# EXAMPLE 1
# The data.frame in data might look like the following:
data<-data.frame(var=paste("X",1:6,sep=""),low=c(1,1,1,0,0,0),
high=c(3,3,3,1,1,1),center=c(2,2,2,0,0,0),nLevels=3,
round=1,factor=0,mix=c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE))
data

# and the design:

optMonteCarlo(~(X1+X2+X3)^2+X4+X5+X6,data)

# Example 2
# Standard designs will often be produced, just as 
# they will with optFederov(). For example,
# a half fraction of a 2^4:
data<-data.frame(paste("X",1:4,sep=""),-1,1,0,2,0,0)
data
optMonteCarlo(~.,data,nTrials=8)

# Example 3
# optMonteCarlo() can treat much larger problems than can 
# optFederov().  For example, optFederov()
# requires a candidate list of 3^20 points for
# a 20 variable, 3 level candidate list -- about
# 25 gigabytes. If the model is quadratic, this must
# be multiplied by about 12. There are other storage
# requirements internal to optFederov() which easily
# double this value. optMonteCarlo() since it only samples
# from the putative candidate list, has no difficulty 
# with a problem of this size. The criterion values
# appearing in the output of optMonteCarlo() are based on
# these samples, but their values seem to be reasonable
# correct, as the following shows: (These are commented
# out for those who have a slow machine.)

dat<-gen.factorial(levels=3,nVar=8)
#desF<-optFederov(~quad(.),dat,eval=TRUE)
#desF[1:5]

data<-data.frame(paste("X",1:8,sep=""),-1,1,0,3,0,0)
#desH<-optMonteCarlo(~quad(.),data,Rand=FALSE,eval=TRUE)
#desH[1:5]

# The following is a 20 variable quadratic. Uncomment
# and wait a while, even if you have a fast machine.
# Note: nRepeats has been changed from its default.
# Note: criterion values for exact designs are often
# far from approximate theory optima; hence, Ge and De
# will be small.

data<-data.frame(paste("X",1:20,sep=""),-1,1,0,3,0,0)
#desBig<-optMonteCarlo(~quad(.),data,nRepeats=1)

# The following will produce improved criterion values

#desNBig<-optMonteCarlo(~quad(.),data,Rand=FALSE,nRepeats=1)

# EXAMPLE 4
# Practically infeasible combinations of variable are 
# common. Designs may be produced which avoid such
# combinations by using a constraint function. Suppose,
# for example that one corner of a cubic box is not
# feasible, then the following will produce a design
# that makes no use of this corner.

Constraints<-function(x){!(x[1]>0.75 && x[2]>0.75)}
data<-data.frame(paste("X",1:4,sep=""),-1,1,0,3,0,0)
desC<-optMonteCarlo(~.,data,con=Constraints)

# The above just removes a corner. Increasing the
# number of levels will remove points along the
# boundary.

data<-data.frame(paste("X",1:4,sep=""),-1,1,0,11,3,0)
desC2<-optMonteCarlo(~.,data,con=Constraints)

}
\keyword{design}
