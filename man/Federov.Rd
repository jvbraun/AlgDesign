\name{optFederov}
\alias{optFederov}
\title{Optimal design}

\usage{
optFederov(frml,data,nTrials,center=FALSE,approximate=FALSE,criterion="D",
	evaluateI=FALSE,space=NULL,augment=FALSE,rows,nullify=0,
	maxIteration=100,nRepeats=5,DFrac=1,CFrac=1,args=FALSE) 

}
\arguments{

	\item{frml}{This may be omitted if \code{data} is the fully model expanded candidate list.
	 If present it should be a formula starting with ~ which describes the model using
	 variable names from \code{data}. It may be ~. if all variables from \code{data} are to be used linearly.
       In addition to the usual operators, quad(), cubic() and cubicS() may be used to 
       expand variables from \code{data} into polynomial models.}
	\item{data}{The candidate list. A matrix or data.frame describing the variables. If a matrix is input and 
	 the columns are not named, they will be assigned names X1,X2, etc. If a data.frame is input without 
       column names, they will be named Var1, Var2, etc. Although data may be input as global variables used 
       in frml, it is preferable to input it here.}
	\item{nTrials}{If approximate=FALSE, it is the number of trials in the final design and if missing, it will be 
	 taken as the greater of \code{length(rows)} or 5 plus the number of terms in the model. If approximate=TRUE,
	 nTrials will be used to round the optimal proportions so that the replications of the points add to nTrials.}
	\item{approximate}{When FALSE, an exact design in nTrails will be calculated. When TRUE the proportions for
	 an approximate theory design will be calculated. If nTrials is set, any proportion less than 1/(2*maxIteration)
	 will be discarded before the proportions are efficiently rounded, otherwise all non-zero proportions will be 	 	 
	 shown: these are the support points.}
	\item{center}{When TRUE, the numeric variables will be centered.}
	\item{criterion}{"D", "A", or "I"}
	\item{evaluateI}{TRUE to evaluate and report I in addition to
	  other criteria. This parameter is included, because evaluating
	  I requires extra effort.}
	\item{space}{If the criterion is "I" or evaluate I is true, the space over which the I criterion is to
	  be evaluated may be input. It should be a matrix with the same column types and names as in data. If
	  space is not input the evaluation will be done over the space described by data.}
	\item{augment}{If TRUE, the row numbers in \code{rows} will never be exchanged.}
	\item{rows}{Either a vector of row numbers (not row names) from \code{data} to be used as the
	  starting design or a vector of row numbers for the design to be augmented. Note, replicate
	  row numbers will be discarded and the length of rows cannot exceed the number of rows in data.}
	\item{nullify}{When non-zero, the initial design is obtained by nullification. If \code{nullify=1},
 	 \code{nTrials} will be calculated (In this case nRepeats is set to 1). If \code{nullify=2}, 
	  number-of-terms trials will be calculated, and the remainder, up to \code{nTrials}, will be filled 
	out at random.}
	\item{maxIteration}{maximum number of times points are exchanged, within each repeat, in seeking an optimum design.}
	\item{nRepeats}{Number of times the entire process is repeated. Has no effect when approximate=TRUE, 
	 or when nullify=1.}
	\item{DFrac}{Design fraction: the fraction of design used in search: 1 uses all of
	  them, 0 uses only the one with the smallest variance.}
	\item{CFrac}{Candidate fraction: the fraction of candidate set searched : 1 uses all of
	  them, 0 uses only the one with the largest variance.}
	\item{args}{If TRUE, the actual arguments to the function including the starting random number
		 seed will be output.}
      }

\description{Calculates an exact or approximate algorithmic design for one of three criteria, using Federov's exchange algorithm.
      }
      
\details{
Let \eqn{E(y)=Zb}, where \eqn{y} is a vector of n observations, \eqn{Z} is an \eqn{n\times k}{n x k} matrix, and \eqn{b} is a vector of k parameters. The ``exact'' design problem is to find a matrix \eqn{Z}, with rows selected from a \eqn{N \times k}{N x k} matrix \eqn{X}, that is ``best'' in some sense. The matrix \eqn{X} can be a discretization of a continuous space or it can represent categories. In either case, the algorithmic design calculation is with respect to \eqn{X}, and not to some larger space containing the points.

Approximate designs weight the candidate points with a probability measure, which for practical purposes amounts to allowing unequal replication. 
  
The Federov(1972) algorithm starts with \eqn{n} points chosen from \eqn{X}. They may be chosen randomly or by nullification, a procedure which iteratively adds points from the null space of \eqn{X}, until a non-singular \eqn{n} point design is found. The Federov algorithm exchanges points in the \eqn{n} point design \eqn{Z} with points in \eqn{X-Z}, i.e. points not in \eqn{Z}, in order to optimize a criterion, and quits when no profitable exchanges are possible, or the input parameter \code{maxIteration} is reached. The quality of the result depends on the starting design and the result may represent a local optimum. The procedure is repeated \code{nRepeats} times in order to come nearer to a global optimum. The parameters \code{DFrac} and \code{CFrac} control the portions of \eqn{Z} and \eqn{X-Z} that are used.


The goal of algorithmic design is to maximize the information about the parameters. The information matrix is a matrix proportional to \eqn{M=Z'Z/n}, and various functions of \eqn{M} are chosen for optimization. The most popular of these is the \code{D} criterion, \eqn{|M|^{1/k}}{|M|^(1/k)}, which is thus a scaling of the ``generalized variance.'' Other criteria of interest involve the variance of predicted values, such as the \code{G} criterion, which is the minimax value of \eqn{d(x)=x'(Mi)x}, over \eqn{X}, where \eqn{Mi} is the inverse of \eqn{M}, and \eqn{x'} is a row of \eqn{X}; and the \code{I} criterion, which is the average value of \eqn{d(x)} in the experimental region. These criteria are invariant under linear transformations of the parameter vector, which frees them from a dependency on units of scale. Other criteria are not invariant, such as the largest eigenvalue of \eqn{M} or the \code{A} criterion, which is \eqn{trace(Mi)/k}: it is of course proportional to the average variance of the parameter estimates. The criteria \code{D}, \code{A}, and \code{I} are supported by \code{optFederov()}, and \code{G}, which is intimately connected to \code{D}, is reported.

The theoretical optimum value of \code{G} is known for approximate theory designs, and so \eqn{G_e}{Ge}, the \code{G} efficiency of \code{G} is available as a standard of design quality. It is especially useful, because \eqn{G_e}{Ge} provides a lower bound on \eqn{D_e}{De}, the \code{D} efficiency for approximate theory, to wit: \deqn{D_e\ge exp(1-1/G_e)}{De >= exp(1-1/Ge)}.

	A vignette giving further details is availble. To access it, type
	
	vignette("AlgDesign")
}

\note{
Algorithmic design is often used with continuous and mixture variables for which R has minimal support, thus the functions quad(), cubic(), and cubicS() may be used in frml. The translation is done with \link{expand.formula}.

Mixture variables are variables such that the rows of \code{data} sum to a constant value, usually unity. Because of the mixture constraint, models containing mixture variables should not have a constant term. The linear and quadratic models for mixture variables A, B, and C are given by \eqn{-1+(A+B+C)} and \eqn{-1+(A+B+C)^2} respectively. See Gorman and Hinman [1962] for 
details.

The function gen.mixture() generates a list of candidate points whose rows sum to unity.


}

\references{
  Atkinson, A.C. and Donev, A.N. (1992). \emph{Optimum experimental
    designs}. Clarendon Press, Oxford.

  Gorman, J.W. and Hinman, J.E. (1962). Simplex lattice designs for 
   multicomponent systems. \emph{Technometrics}. 4-4. 463-487.

  Federov, V.V. (1972). \emph{Theory of optimal experiments}. Academic
   Press, N.Y.
  

}

\value{

  \item{D}{The kth root of the generalized variance: \eqn{det(M)^{1/k}}{det(M)^(1/k)}, where \eqn{det(M)} is the
    determinant of the normalized dispersion matrix \eqn{M} --
    i.e. \eqn{M=Z'Z/n}, where \eqn{Z=X[rows,]}}
  \item{A}{The average coefficient variance: \eqn{trace(Mi)/k}, where \eqn{Mi} is the inverse of \eqn{M}.}
  \item{I}{The average prediction variance over X, which can be shown to be \eqn{trace((X'X*Mi)/N)}, where
	N is the number of rows in X. This is calculated only when I is the criterion or when \code{evaluateI} is TRUE.}
  \item{Ge}{The minimax normalized variance over X, expressed as an efficiency with respect to the optimal approximate
	theory design. It is defined as \eqn{k/max(d)}, where \eqn{max(d)} is the maximum normalized
    	variance over \eqn{X} -- i.e. the max of \eqn{x'(Mi)x}, over all rows \eqn{x'}
    	of \eqn{X}.}
  \item{Dea}{A lower bound on \code{D} efficiency for approximate theory designs. It is equal to \eqn{exp(1-1/Ge)}.}
  \item{design}{The design.} 
  \item{rows}{A numerical vector of the design row numbers.}
  \item{args}{A list of the actual arguments used in this call.}
}


\examples{

# EXAMPLE 1
# A quadratic polynomial in three variables. The resulting D will be about 0.46.
# This may be compared with a standard central composite design obtained from
# rows 1,3,5,7,9,11,13,15,17,19,21,23,25,27 of dat, which has a D value of 0.46.
# The central composite design seems to be the optimal design for all three criteria.

dat<-gen.factorial(levels=3,nVars=3,varNames=c("A","B","C"))

desD<-optFederov(~quad(A,B,C),dat,nTrials=14,eval=TRUE)
desA<-optFederov(~quad(.),dat,nTrials=14,eval=TRUE,crit="A")
desI<-optFederov(~quad(.),dat,nTrials=14,eval=TRUE,crit="I")

rows<-c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)
desO<-optFederov(~quad(.),dat,nTrials=14,eval=TRUE,rows=rows)

# The I criterion may be seen to decrease as the space is expanded. 

levels<-seq(-1,1,by=.1)
dat<-expand.grid(list(A=levels,B=levels,C=levels))

desL<-optFederov(~quad(.),dat,nTrials=14,eval=TRUE)

# This is not the case for A or D. For A and D, the support points are the points 
# of the grid with the three levels above. Points not on this grid move
# the criteria in a non-optimal direction; hence, the enlarging space has no effect.

# EXAMPLES 2
# Standard designs are usually optimal designs. If nTrials is set to that for
# a standard design, and if nRepeats is large enough, the standard design will 
# often be found For example, a half replicate of a 2^4 will be obtained by the 
# following. 

dat<-gen.factorial(levels=2,nVars=3,varNames=c("A","B","C"))
desH<-optFederov(~.,dat,8)

# A third replicate of a 3^3 will be obtained by the following:

dat<-gen.factorial(levels=3,nVars=3,factor=1:3)
desT<-optFederov(~.,dat,9)

# An orthogonal design similar to a 12 run Plackett-Burman design can be 
# created by the following. 

dat<-gen.factorial(levels=2,nVars=11,varNames=c("A","B","C","D","E","F","G","H","J","K","L"))
desPB<-optFederov(~.,dat,12,nRepeats=20)

# The above calculation is numerically difficult for the A and I criteria, 
# and nRepeats=100 or more may be needed. 

# It is instructive to examine a case in which the standard design is not found.
# The following is an attempt to create a Latin square design. It is not always successful.

lv<-factor(1:5)
dat<-expand.grid(A=lv,B=lv,C=lv)
desL<-optFederov(~.,dat,nTrials=25,nRep=100)

# It may be summarized as follows.

cs<-xtabs(~.,desL$design)
{xx<-matrix(0,5,5); for (i in 1:5) xx=xx+cs[1:5,1:5,i]*i;xx}
 


# EXAMPLE 3
# Mixture variables have a constant sum, usually 1. This causes a linear dependency
# among terms of polynomial models. In particular the constant term is dependent.
# Squared terms in a quadratic model are confounded with interaction terms, so that
# a quadratic model for three mixture variables is ~0+(A+B+C)^2. The following
# calculation generates a set of candidate varibles using gen.mixture() with
# four values on each axis, and then creates a 15 run design. The design is optimal.
# Indeed, the candidate set produced by gen.mixture(2,5) is optimal. Note: 
# nullify=TRUE is used to ensure that this example will run withough error. The
# default value of 5 for nRepeats is sometimes not enought to find a starting
# design with a mixture problem.


dat<-gen.mixture(4,5)
desM<-optFederov(~(X1+X2+X3+X4+X5)^2-1,dat,15,nullify=TRUE)

# EXAMPLES 4
# Design augmenation can be obtained by setting augment=TRUE, and placing the row numbers
# of the design to be agmented in rows. Augmentation is often used to (1) add a new variable
# to an existing design or (2) to increase the complexity of the model. The following illustrates
# adding a variable to an existing design using desD above. It is assumed that all runs of the
# existing design have been made at the -1 level of the new variable:

dat<-gen.factorial(levels=3,nVars=3,varNames=c("A","B","C"))
desA<-optFederov(~quad(.),dat,nTrials=25,augment=TRUE,rows=desD$rows)

# The half fraction in desH, can be augmented to support an additional term:

dat<-gen.factorial(levels=2,nVars=4,varNames=c("A","B","C","D"))
desH<-optFederov(~.,dat,8)
desH2<-optFederov(~A+B+C+D+I(A*B),dat,10,aug=TRUE,rows=desH$rows)

# EXAMPLES 5
# Optimal approximate theory designs have non-zero probabilities only on support points.
# For the first example above the approximate theory design is as follows. It shows
# that all points in the cubic lattice are support points. The D for this 
# design is 0.474 which may be compared with the D of 0.463 of the first example to
# indicate that that exact design had a D-efficiency of 97\%. The lower bound Dea
# was 82\%.

dat<-gen.factorial(levels=3,nVars=3,varNames=c("A","B","C")) 
desDA<-optFederov(~quad(A,B,C),dat,eval=TRUE,approx=TRUE)

# The largest proportions will be rounded if nTrials is specified.

desDAN<-optFederov(~quad(A,B,C),dat,eval=TRUE,approx=TRUE,nTrials=15)
}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
	
	Please cite this program as follows:

Wheeler, R.E. (2004). optFederov. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}


\keyword{design}
