\name{optBlock}
\alias{optBlock}
\title{Optimal design blocking}

\usage{

optBlock(frml,withinData,blocksizes,rows=NULL,wholeBlockData=NULL,center=FALSE,
	nRepeats=5,criterion="D",args=FALSE)

}
\arguments{
	\item{frml}{This may be omitted if \code{data} is the fully model expanded candidate list.
	 If present it should be a formula starting with \code{~} which describes the model using
	 variable names from \code{data}. It may be \code{~}. if all variables from \code{data} are to be 
	 used linearly. It may also be \eqn{~.^p}, where p is an integer. In addition to the usual 
	 operators, quad(), cubic() and cubicS() may be used to expand variables from \code{data} 
	 into polynomial models.}
	\item{withinData}{A matrix or data.frame describing the variables. If the columns are not 
       named, they will be assumed to have the names X1,X2, etc. Although data may be input as
	 global variables used in frml, it is preferable to input it here. The number of rows in \code{withinData}
	 must be at least as large as the sum of the number of terms plus the number of blocks.}
	\item{blocksizes}{A vector giving the block sizes for each block. The length of \code{blocksizes}
		specifies the number of blocks.}
	\item{rows}{Row numbers (not rownames) of withinData rows to be used as a starting blocked design.}
	\item{wholeBlockData}{A matrix or data.frame describing the whole block variables. The number of
		rows should equal the length of blocksizes. Each row should
		correspond to the settings of the whole block variable in a block.}
	\item{center}{If TRUE, the withinData and wholeBlockData will be centered.}
	\item{nRepeats}{The number of times the entire process is repeated.}
	\item{criterion}{"D":D-criterion: "OB","OBS": orthogonal blocks, unscaled and scaled; "Dp","Dpc": 
		Product of block determinants, uncentered and centered.}
	\item{args}{If TRUE, the actual arguments to the function including the starting random number
		 seed should be output.}
      }

\description{Blocking of experimental designs using various criteria.}
      
\details{
Cook and Nachtsheim (1989) developed a blocking algorithm using the D criterion, which not only finds blocks, but attempts to optimize the design itself. The package Dopt, converted by V.N. Venables, implements this technology -- unfortunately there seem to be bugs in the FORTRAN code. 

This function uses various criteria to block either pre-existing designs or to creates a new blocked design from a candidate list. There may be interactions between whole plot and within plot factors, and it may be used to produce multistratum blocked designs. 

The algebraic expressions used here for the D criterion are quite different than those used by either Cook and Nachtshim or Goos and Vandebroek. They are described in \emph{Comments on algorithmic design}, a paper accompanying this package.

Although the D criterion produces good designs for all blocks, the allocation for individual blocks can often be improved upon by the use of either the Dp or Dpc criterion, which optimize the product of the determinants of the individual blocks. The Dp criterion uses uncentered blocks, the Dpc uses centered blocks. 

Blocking may also be done using the orthogonal blocking procedure of Nguyen. This comes in two flavors OB and OBS which differ in that for OBS, the columns of S are scaled by division with the sample variance of the data, thus deemphasizing squared and other large terms.

Blocking may be done when whole plot factors interact with within plot factors. Split plot experiments are sometimes of this nature. Goos and Vandebroek (2003) developed a procedure for this problem which algorithmically blocks using a candidate list and the ratio between the whole and within variances. Trinca and Gilmour (2001) give an algorithm which does not depend on the ratio of the variances. The present procedure assumes that the whole blocks and their factors, have been decided upon in advance. It makes no assumptions about the ratio between whole and within variances. 

	A vignette giving further details is availble. To access it, type
	
	vignette("AlgDesign")

}



\references{
Atkinson, A.C. and Donev, A.N. (1989). The construction of exact D-optimum experimental designs with application to blocking response surface designs. \emph{Biometrika}. 76. 515-526.

Cook, R.D. and Nachtsheim, C.J. (1989). Computer-aided blocking of factorial and response-surface designs. \emph{Technometrics.} 31-3. 339-346.

Goos, P. and Vandebroek, M. (2003). D-optimal split-plot designs with given numbers and sizes of whole plots. \emph{Technometrics.} 45-3. 235-245.

Nguyen, Nam-Ky. (2001). Cutting experimental designs into blocks. \emph{AusNZJSt}. 43-3. 367-374.

Trinca, L.A. and Gilmour, S.G. (2000). An algorithm for arranging response surface designs in small blocks. \emph{Computational Statistics and Data Analysis}. 33. 25-43.

Trinca, L.A. and Gilmour, S.G. (2001). Multistratum response surface designs. \emph{Technometrics}. 43-1. 25-33.

}

\examples{
# Blocking the design for a quadratic polynomial in three variables into two 
# seven trial blocks:

dat<-gen.factorial(3,3,varNames=c("A","B","C"))
desD<-optFederov(~quad(.),dat,nTrials=14,eval=TRUE) # Choose an optimum 14 trail design.
optBlock(~quad(.),desD$design,c(7,7))

# Letting optBlock() search the dat candidate list instead of first choosing a 
# 14 trial design.
optBlock(~quad(.),dat,c(7,7))


# A block design for 7 treatments in 7 blocks of size 3. Note how withinData 
# is recycled to fill out the blocksize requirements.

BIB<-optBlock(~.,withinData=factor(1:7),blocksizes=rep(3,7))

# This is a balanced incomplete block design as may be seen from:

crossprod(table(BIB$rows,c(rep(1:7, rep(3,7)))))

# A partially balanced incomplete block design with two associate classes:

tr<-factor(1:9)
PBIB<-optBlock(~.,withinData=tr,blocksizes=rep(3,9))

crossprod(table(PBIB$rows,c(rep(1:9, rep(3,9)))))


# Two fractions of a 2^(4-1).

dat<-gen.factorial(2,4)
od<-optBlock(~.,dat,c(8,8)) 

# The blocks are not themselves orthogonal even though the entire design is optimal.

bk<-data.matrix(od$Blocks$B1)
t(bk)\%*\%bk 

# Better blocks may be obtained as follows, but note that they are not generally 
# the fractions that would be obtained by confounding the third order interaction.

od<-optBlock(~.,dat,c(8,8),criterion="Dpc",nR=10)
bk<-data.matrix(od$Blocks$B1)
t(bk)\%*\%bk

# Blocking with whole plot factors. Note that the 27 rows of within are recycled 
# to make the 54 trial blocked design. 

within<-expand.grid(A=c(-1,0,1),B=c(-1,0,1),C=c(-1,0,1))
whole<-expand.grid(D=factor(1:3),E=factor(1:3))
od<-optBlock(~D+E*(quad(A,B,C)),withinData=within,blocksizes=rep(6,9),wholeBlockData=whole)

# Either withinData, or wholeBlockData may be an approximate theory optimial design 
# produced by optFederov() for nTrials. The first column in the optFederov() output 
# design, named "Rep..", is used to replicate the trials. 

within<-optFederov(~quad(A,B,C),within,nT=54,approx=TRUE)
od<-optBlock(~D+E*(quad(A,B,C)),withinData=within$design,blocksizes=rep(6,9),wholeBlockData=whole)

}


\value{

  \item{D}{\eqn{det(M)^{1/k}}{det(M)^(1/k)}, where \eqn{det(M)} is the
    determinant of the normalized dispersion matrix \eqn{M} --
    i.e. \eqn{M=X'X/N}, where each row of X has had the appropriate block mean subtracted.}
  \item{Dp (Dpc)}{\eqn{\left(\prod_1^b (det(M_i))^{1/k}\right)^{1/b}}{[prod (det(Mi)^(1/k)]^(1/b)}, 
	where \eqn{det(M_i)}{det(Mi)} is the determinant of the normalized dispersion matrix 
	for the uncentered (centered) block submatrix of \eqn{X}, and \eqn{b} is the number of blocks. The
	normalization is with respect to the number of observations in each block.}
  \item{diagonality or SS}{The diagonality is \eqn{(|M|/P)^{1/k}}, where P is the product of the diagonal
	elements of M. The SS is the sum of squares of S when the OB criterion is used.}
  \item{Blocks}{A list of the blocks, labeled B1, B2, etc.}
  \item{design}{A data.frame. The design with blocks stacked in order.}
  \item{rows}{Numeric row numbers of the design rows corresponding to the withinData rows.} 
  \item{args}{A list of the actual arguments used in this call.}

}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}

	Please cite this program as follows:

Wheeler, R.E. (2004). optBlock. \emph{AlgDesign}. The R project for statistical computing \url{https://www.r-project.org/}
}
\keyword{design}

