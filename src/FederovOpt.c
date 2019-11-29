/* file FederovOpt.c
| copyright (C) 2002-2004 by Robert E. Wheeler
*/


#include "wheeler.h"

#include <string.h>
#include <stdlib.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Utils.h>


bool DllMain(void)					 
{

    return true;
}

/***********************************************************************/

#define Imat(R,C)  ((C)+((R)*(2*nColumns-(R)-1))/2)  /* Upper triangular plus rect block*/

int		Klimit=0;
int		Llimit=0;
bool    doRcrit=FALSE;

typedef struct tagdStruct{
	double	d;
	int		i;	/* item index */
	int		o;	/* d sorted index -- smallest d has index 0 */
} dType, *pdType;

bool doApprox=false; /* global switch */

int FederovOptimize(double *X,double *B,double *BU,double *proportions,bool RandomStart,int Nullify,
	int	criterion,bool evaluateI,bool doSpace,
	int augment,double *D,double *A,double *I,double *G,double *U,double *V,double *T,double *Ti,double *Tip,
	double *W,double *maxmin,dType *d,double *vec,int *designFlag,int *ttrows,int *rows,
		    int	*trows,int N,int n,int k,int maxIteration,int nRepeats,double DFrac,double  CFrac,bool *error);


SEXP FederovOpt(SEXP Xi,SEXP RandomStarti,SEXP rowsi,SEXP Nullifyi,SEXP criterioni,SEXP evaluateIi,
		SEXP doSpacei, SEXP Bi,SEXP augmenti,SEXP doApproxi,SEXP proportionsi,SEXP ni,SEXP maxIterationi,SEXP nRepeatsi,SEXP DFraci,SEXP CFraci);

double GetLinearCriterion(double *pBi,double *pBj,int criterion,double *pVi,double *pVj,
	double di,double dj,double dij,int k,double dn);
double GetLinearCriterionA(double *pBU,int criterion,double *pU,int k);

int  ProgAlloc(double **U,double **V,double **B,
	double **BU,double **T,double **Ti,double	**Tip,double  **W,double **maxmin,dType **d,
	double  **vec,int **designFlag,int **ttrows,int **trows,int Nin,int	n,int k,bool criterion,
	bool evaluateI,bool doSpace);
void ProgDealloc(double *U,double *V,double *B,double *BU,double	*T,double *Ti,
	double	*Tip,double  *W,double *maxmin,dType *d,double  *vec,bool *designFlag,bool *ttrows,int *trows);
void FillInB(double *X,double *B,int k,int N);

void BacksolveT(double *matrixXY,int nColumns,bool scaled);
double findDelta(double *B,int criterion,int *xold,int *xnew,dType *d,double *U,double *V,double N,int n,int k,
	bool *designFlag,int *rows,bool *failure);
double findDeltaAlpha(double *delta,double *BU,int criterion,int *xnew,double maxd,int maxdi,
	double Acrit,double Icrit,dType *d,
	double *U,double N,int k,bool *failure);
void getEfficiencies(double	*determinant,double	*Gefficiency,double	*DRefficiency,
	bool *designFlag,double	*X,double *T,double	*vec,dType	*d,int m,int N,int k);
double  makeTiAndTipFromT(double *Tip,double *T,double *Ti,double *maxmin,double norm,bool *singular,int k);
double evaluateCriteria(double *Tip,double *Ti,double *W,double *B,int criterion,bool evaluateI,
	double  *Acrit,double  *Icrit,double logdet,int k,int N);
void makeBXd(double *X,double *U,double *V,double *Ti,double *Tip,double *B, double *BU, int criterion,
			 int *designFlag,dType *d,double *maxd,int *maxdi,int k,int N);
void MatMult(double *A,double *B,double *C,int a,int b);
void Permute(int *a,int	n);
void getRange(double *pMx,double *pMn,double *pX,int k);
double reduceDesign(int	*rows,double *X,double *T,double *maxmin,double *tVec,int k,
		int n,bool scale,bool *singular);
void Rotate(double	*vec,double	*tVec,double *matrixXY,int nTerms,int nColumns,double weight,double n);
void SwapRows(pdType a,pdType b);
int  dCompare(pdType a,pdType b,bool compType);
void dShellSort(pdType pd,int n,bool compType);
void update(int	xold,int xnew,int *rows,int *designFlag,double *T,double *X,double *tVec,
	int k,int n);
void updateA(int xnew,double *proportions,double alpha,double *T,double	*X,double *tVec,int	k,int N);
double getNextRow(double *V,int N,int k,int *designFlag,int	*newRow);
void orthog(double *V,double *vec,int *designFlag,double scale,int N,int k);
void orthogAug(double *V,int *rows,int augment,int *designFlag,int N,int k);
int nullify(double *X,double *V,int augment,int *rows,bool *designFlag, int N,int k);
void filloutDesign(double	*T,double *Ti,double *Tip,int	*rows,
	int *ttrows, double	*X,double *maxmin,double *vec,
	int k,int ka,int	n,int N,bool *singular);
void transposeMatrix(double *X,int N,int k);


/* dShellSort
|  Sorts an array of pdType in increasing order
*/


void SwapRows(
	pdType	a,
	pdType	b
)
{
   dType temp;

   temp=*a;
   *a=*b;
   *b=temp;
}

int dCompare(	
	pdType	a,
	pdType	b,
	bool		compType /* true for d, false for i */
)
{
	double	d;
	int		c;

	if (compType) {
		d=a->d-b->d;
		if (d<0)
			return(-1);
		if (d==0)
			return(0);
		else
			return(1);
	}
	else {
		c=a->i-b->i;
		if (c<0)
			return(-1);
		if (c==0)
			return(0);
		else
			return(1);

	}
}

void dShellSort(
	pdType	pd,		/* the array of dType */
	int		n,       /* size of the array     */
	bool	compType /* true if d's compared, false if i's   */
)
{
   int		gap,
			i,
			j;
   pdType	p,
			q;

   for(gap=n/2;gap>0;gap/=2)
   {
		for(i=gap;i<n;i++)
		{
			for(j=i-gap;j>=0;j-=gap)
			{
				p=pd+j;
				q=pd+j+gap;
				if (dCompare(p,q,compType)<=0)
					break;
				SwapRows(p,q);
         }
		}
	}
}




/* Permute **********************************************************
|	Randomly pemutes the n integers in a[] using the Fike
|	algorithm.  See Fike, "A permutation generation method"  The Computer
|	Journal, 18-1, Feb 75, 21-22.
*/

void Permute(
int		*a, 
int		n
)
{
   int 		i,
			j,
			temp;

   GetRNGstate();
   for (i = 1; i < n; i++) {
      j =(int)((double)(1+i)*unif_rand());
      temp = a[j];
      a[j] = a[i];
      a[i] = temp;
   }
   PutRNGstate();
}


/* Rotate ****************************************************************
|  Rotates in the new row using Given's rotations without
|  square roots.  weight is the scale factor for the new row.
|  For X'(X,Y)=T'D(T,U), the result (returned in matrixXY) is (T,U) with unit
|   diagonals.  These being understood, D is placed on the diagonal.
|	NOTE: matrixXY is nTerms by nColumns, with the T part stored as an upper triangular matrix
|		requiring nTerms(nTerms+1)/2 elements. Imat indexes the elements.
|	tVec is needed if vec is to be preserved.
|	If weight is -1, vec is removed from T
|	NOTE: the effect is to multiply each row by sqrt(weight)
|	matrixXY should initially be zero unless vec is to be added or removed from existing T
*/

#define TOLROT 1.0e-50

void Rotate(
	double	*vec,
	double	*tVec,
	double	*matrixXY,
	int		nTerms,
	int		nColumns,
	double	weight,
	double  n
)
{
	double	d,
			dn=sqrt(n),
			dp,
			c,
			s,
			x,
			r;
	int		i,
			j,
			kIndex;
	bool	skip;


	for (i=0;i<nColumns;i++) {
		tVec[i]=vec[i]/dn;       /* dividing by sqrt(n) to produce a normalized result */
								 
	}

	skip=false;
	for (i=0;i<nTerms;i++) {
		if(! skip) {
			x=tVec[i];
			if (x==0.0)
				continue;
				
			d=matrixXY[kIndex=Imat(i,i)];
			dp=d+weight*x*x;
			if (fabs(dp)<TOLROT)
				continue; 
			matrixXY[kIndex]=dp;
			c=d/dp;
			s=weight*x/dp;
			
			if(d==0.0) {
				skip=true;
				weight=0.0;
			} 				/* to avoid 0/0 since matrixXY is originally 0*/
			else 
				weight*=c;				
			kIndex++;
			for (j=i+1;j<nColumns;j++,kIndex++)	{
				r=matrixXY[kIndex];
				matrixXY[kIndex]=s*tVec[j]+c*r;
				tVec[j]-=x*r;
			}
		}
		else
		  break;
	}


}

/* BacksolveT **************************************************************
|  Backsolves T in place to obtain T inverse. T is stored in matrixXY with T
|		upper triangular.  See Rotate() for a description of T.
|  The dimensions of is nColumns by nColumns
|	If scaled is true, the rows are scaled. If scaled is false, the diagonal
|	contains the scale values and the diagonals of T are understood to be unities.
*/

void BacksolveT(
	double	*matrixXY,
    int		nColumns,
	int		scaled
)
{
   int	col,
	    i,
		j,
		lIndex,
		kIndex;


	if (scaled) {
 		for (col=nColumns-1;col;col--) {
			matrixXY[Imat(col,col)]=1/matrixXY[Imat(col,col)];
			for (j=col;j--;) {
			   lIndex=Imat(j,col);
			   kIndex=Imat(j,j+1);
				matrixXY[lIndex]*=-matrixXY[Imat(col,col)];
				for (i=1;i<col-j;i++) {
					matrixXY[lIndex]-=matrixXY[Imat(i+j,col)]*
						matrixXY[kIndex++];
				}
				matrixXY[lIndex]/=matrixXY[Imat(j,j)];
			}
		}
		matrixXY[0]=1/matrixXY[0];
	}
	else {
			/* Backsolve for the inverse of T, which is upper triangular */
 		for (col=nColumns-1;col;col--) {
			for (j=col;j--;) {
			   lIndex=Imat(j,col);
			   kIndex=Imat(j,j+1);
				matrixXY[lIndex]*=(-1.0);
				for (i=1;i<col-j;i++)
					matrixXY[lIndex]-=matrixXY[Imat(i+j,col)]*
						matrixXY[kIndex++];
			}
		}

		for (i=0;i<nColumns;i++)
			matrixXY[Imat(i,i)]=1.0/matrixXY[Imat(i,i)];
	}

}




/* getRange ***********************************************************************
|  increments max and min with values from a new vector
*/

void getRange(
	double *pMx,
	double *pMn,
	double *pX,
	int k
)
{
	int i;

	for (i=0;i<k;i++) {
		pMx[i]=maxm(pMx[i],pX[i]);
		pMn[i]=minm(pMn[i],pX[i]);
	}

}




/* reduceDesign ********************************************************************
|	Reduces design rows in X to T via Given's rotations. Returns the determinant
*/

double reduceDesign(
	int		*rows,
	double	*X,
	double	*T,
	double  *maxmin,
	double	*tVec,
	int		k,
	int		n,
	bool    scale, /* if true, scale T by n */
	bool    *singular
)
{
	double	*pX;
	double  *pMx=maxmin;
	double  *pMn=maxmin+k;
	double	logdet=0;
	double  maxT=-1e16;
	double  minT=1e16;
	double  t;
	double  r;
	double  norm;
	int		K=k*(k+1)/2;
	int		sizeT=K*sizeof(double);
	int		nColumns=k; /* Needed by Imat() */
	int		i;

	*singular=false;

	for (i=0;i<k;i++) {
		pMx[i]=maxT;
		pMn[i]=minT;
	}

	memset((void *)T,0,sizeT);


	norm=scale?(double)n:1.0;
	for (i=0;i<n;i++) {
		pX=X+rows[i]*k;
		getRange(pMx,pMn,pX,k);
		Rotate(pX,tVec,T,k,k,1.0,norm);  /* the n normalizes the design */
	}

	norm=2*scale?sqrt((double)n):1.0;
	for (i=0;i<k;i++) {
		r=(pMx[i]+pMn[i])/norm;
		t=T[Imat(i,i)];
		if (t<=0  || t<r*1e-16) {
			*singular=true;
			return(0);
		}
		logdet+=log(t); 
	}


	return (logdet);
}

/* getNextRow ********************************************************************
|  Returns index of the longest row of V in *newRow and its squared length
*/

double getNextRow(
	double  *V,
	int		N,
	int		k,
	int		*designFlag,
	int		*newRow
	)
{	
	int i;
	int j;
	double scale=-1;
	double ssq;
	double *pV;

	*newRow=-1;
	for (i=0;i<N;i++) {
		if (!designFlag[i]) {
			ssq=0;
			pV=V+k*i;
			j=k;
			while(j--) {
				ssq+=(*pV)*(*pV);
				pV++;
			}
			if (ssq>scale){
				*newRow=i;
				scale=ssq;
			}
		}
	}
	

	return(scale);
	
}


#define nullTol 1e-8


/* orthog ***************************************************************************
|	Orthogonalizes V with respect to newRow
*/

void orthog(
	double	*V,
	double  *vec,
	int		*designFlag,
	double  scale,
	int		N,
	int		k
	)
{
	int i;
	int j;
	double beta;
	double cross;
	double *pvec;
	double *pV=0;
	double *ppV;



	for (i=0;i<N;i++) {
		pV=V+i*k;
		if (!designFlag[i]) {
			cross=0;
			pvec=vec;
			j=k;
			ppV=pV;
			while(j--) {
				cross+=(*ppV++)*(*pvec++);
			}
			beta=cross/scale;
			ppV=pV;
			pvec=vec;
			j=k;
			while(j--) {
				(*ppV++)-=beta*(*pvec++);
			}		
		}
	}	
}


/* orthogAug ***************************************************************************
|	orthogonalies V with respect to augment rows
*/

void orthogAug(
	double	*V,
	int		*rows,
	int		augment,
	int		*designFlag,
	int		N,
	int		k
	)
{
	int i;
	int j;
	double *pV;
	double *ppV;
	double scale;

	for (i=0;i<augment;i++) {
		ppV=pV=V+rows[i]*k;
		scale=0;
		j=k;
		while(j--) {
			scale+=(*ppV)*(*ppV);
			ppV++;
		}
		orthog(V,pV,designFlag,scale,N,k);		
	}

}

/* nullify *************************************************************************
|  Nullification algorithm
|	Finds next row to add as longest in null space of existing design
|	 then orthogonalizes space with respect to this new row, and repeats 
|	 till k rows have been found or singular design.
*/

int nullify(
	double	*X,
	double	*V,
	int		augment,
	int		*rows,
	bool    *designFlag,
	int		N,
	int		k
	)
{
	int i;
	int newRow;
	int sizeT=N*k*sizeof(double);
	double scale;
	double tolerance=nullTol;

	memcpy((void*)V,(void*)X,sizeT);

	if (augment) {
		orthogAug(V,rows,augment,designFlag,N,k);
	}

	for (i=0;i<k-augment;i++) { 
		scale=getNextRow(V,N,k,designFlag,&newRow);
		if (scale<tolerance || newRow==-1)
			return(0); /* singular design. */
		if (i==0)
			tolerance=scale*nullTol;

		rows[i+augment]=newRow;
		designFlag[newRow]=true;

		if (i!=(k-1)) {
			orthog(V,V+newRow*k,designFlag,scale,N,k);
		}
	}

	return(1);
}

#define detTol 1e-16

/* makeTiAndTipFromT ******************************************************************
|	Finds Ti, (upper triangular) the inverse of T, and then places Ti' (lower triangular) into Tip
|	NOTE: T has scale values on diagonal which are incorporated in Tip
|	returns logdet
*/

double  makeTiAndTipFromT(
	double *Tip,
	double *T,
	double *Ti,
	double *maxmin,
	double norm,
	bool   *singular,
	int		k
	)
{
	double logdet=0;
	double *pTi;
	double d;
	double r;
	double *pMx=maxmin;
	double *pMn=maxmin+k;
	int		g;
	int K=k*(k+1)/2;
	int	sizeT=sizeof(double)*K;
	int		i;
	int		j;
	int	nColumns=k;	/* needed by Imat(), which indexes the upper triangle of a matrix. */


	*singular=false;

	memcpy((void*)Ti,(void*)T,sizeT);

	/* Scale rows of Ti by diagonal and get log determinant */
	pTi=Ti;
	for (i=0;i<k;i++) {
		r=(pMx[i]-pMn[i]);  /* The design really shouldn't be singular, but in extreme */
		r=detTol*r*r/norm;  /* cases the A and I criteria can fail numerically. */
		if (*pTi>r)		
			logdet+=log(*pTi);
		else
			*singular=true;
		d=sqrt(*pTi);
		(*(pTi++))=d;
		for (j=i+1;j<k;j++)
			(*pTi++)*=d;
	}


	BacksolveT(Ti,k,TRUE); /* Get the inverse of T in Ti */

	/* Put transpose of Ti in Tip */
	g=0;
	for (j=0;j<k;j++) {
		for (i=0;i<=j;i++) {
			Tip[g++]=Ti[Imat(i,j)];
		}
	}

	return(logdet);

}

/* evaluateCriteria **********************************************************************
|	The trace of TiTi'=Tip'Tip is returned as Acrit, and the trace of G*TiTi' in Icrit, where G
|		is a symmentric matrix whose upper triangle, including the diagonal, is in B.
*/

double evaluateCriteria(
	double	*Tip,
	double  *Ti,
	double	*W,	/* k*(k+1)/2 elements of working strorage */
	double  *B, /* Upper trianglular part of symmetric matrix G */
	int		criterion,
	bool	evaluateI,
	double  *Acrit,
	double  *Icrit,
	double  logdet,
	int		k,
	int		N
)
{	int	i;
	int j;
	int g;
	int K=k*(k+1)/2;
	int	sizeT=sizeof(double)*K;
	double t;
	double *pTip;
	double *pTip0;
	double *pW;
	double *pB;
	double *pTi;
	int	nColumns=k;	/* needed by Imat(), which indexes the upper triangle of a matrix. */

	/* Find the trace of TiTi' and put it in *Acrit */
	
	t=0;
	pTi=Ti;
	for (i=0;i<K;i++) {
		t+=(*pTi)*(*pTi);
		pTi++;
	}

	*Acrit=t/(double)k;

	if (criterion==2|| evaluateI) {
		/* To find the I criterion which is trace G*TiTi'. */
		/* Premultiply column i of TiTi' by row i of G, and then sum to get the trace. */
		*Icrit=0;
		for (i=0;i<k;i++) { /* Cycle through the rows of G */
			/* First form the ith column of TiTi' in W. */
			memset((void *)W,0,sizeT);
			pTip0=Tip+(i*(i+1))/2;           /* First element of row i of Ti' */
			pTip=pTip0+i;	/* last element of row i of Ti' == first non-zero element in column i of Ti'*/
			for (j=i;j<k;j++) { /* Successive elements in column I of Ti'*/
				pW=W;   /*Start storage with first element of W */
				for (g=0;g<=j;g++) { /* Cycle through the non-zero elments in row j of Ti' */
					(*pW++)+=(*pTip0++)*(*pTip);
				}
				pTip=pTip0+i; /* Next element in column i of Ti' */
			}
			/* Now premultply by G and sum the products in t */
			t=0;
			for (j=0;j<i;j++) { /* First the i elements in column i above the diagonal */
				t+=B[Imat(j,i)]*W[j];	
			}
			/* Then the remaining elements in row i from the diagional onward */
			pB=B+Imat(i,i);
			for (j=i;j<k;j++) {
				t+=(*pB++)*W[j];
			}
			*Icrit+=t;
		}
	}

	switch (criterion) {
	case 0:
		return (logdet); /* Simply return the current determinant */
		break;
	case 1:
		return (*Acrit);
		break;
	case 2:
		return (*Icrit);
		break;
	}
	return 0;
}

/* MatMult *****************************************************************************
|  Matrix multiply, Tip*B=C
|	A is the upper triangle of a symmetric matrix Tip of dim axa
|	B is axb, but B' of dim bxa is stored.
|	C is axb, but C' of dim bxa is stored.
*/

void MatMult(
	double *A,
	double *B, /* actually B' */
	double *C, /* actually C' */
	int		a,
	int		b
	)
{	
	int i;
	int j;
	int k;
	double *pA;
	double *pB;
	double *pC;
	double *pBj;
	double s;

	memset((void *)C,0,sizeof(double)*a*b); /* Zero C */

		/* Multiply C=AB */
	pB=B;	/* First column of B and C */
	pC=C;
	for (i=0;i<b;i++) { /* Cycle through all columns of B and C */
		pA=A;			/* Reset pA to the first element of A */
		for (j=0;j<a;j++) { /* Cycle through each row of A */
			pBj=pB+j;	/* Get next B element */
			s=0;
			for (k=j;k<a;k++) { /* Multiply  row of A by column of B */
				s+=(*pA++)*(*pBj++);
			}
			(*pC++)=s; /* Store the result in C */
		}
		pB+=a; /* Next column of B */
	}
	/* Multiply C=A1'B+C [=(A1'+A)B], Where A1' is A' without a diagional */
	pB=B; /* First column of B and C */
	pC=C;
	for (i=0;i<b;i++) { /* Cycle through all columns of B and C */
		pA=A;				/* Reset pA to the first element of A */
		pBj=pB;		/* First element of current column of B */
		for (j=0;j<a-1;j++) { /* Cycle through the first a-1 elements of a column of B */
			pA++;		/* Skip the diagonal element of A */
			for (k=j+1;k<a;k++) { /* Multiply the element of B times elements in the row of A1' and  */
								  /* add the product to the corresponding element of C */
				(*(pC+k))+=(*pA++)*(*pBj);  /* Note: elements of C are cycled through */
			}
			pBj++; /* Next element in current column of B */
		}
		pB+=a; /* Next column of B and C */
		pC+=a;
	}
}

/* makeBXd ***************************************************************************
|	Makes V=XTi, U=XTiTi'  and d, which is a vector comprising sums of squares of V rows
|	Then makes BU if criterion==2
*/

void makeBXd(
	double	*X,
	double	*U,
	double  *V,
	double  *Ti,
	double	*Tip,
	double  *B,
	double  *BU,
	int		criterion,
	int		*designFlag,
	dType	*d,
	double  *maxd,
	int		*maxdi,
	int		k,
	int		N
)
{
	int		i;
	int		j;
	int		l;
	double	*pX=X;
	double	*pV;
	double  *ppV;
	double  *pU;
	double	sum;
	double	sumsq;
	double	*p;
	double	*q;
	dType	*pd=d;
	int		g;
	
	for (i=0;i<N;i++) {
		pV=V+i*k;
		q=Tip;
		sumsq=0;
			/* making row i of XTi */
		ppV=pV;
		for (j=0;j<k;j++) {
			p=pX; /* row i of X */
			sum=0;
			for (l=0;l<=j;l++) {
				sum+=(*p++)*(*q++);	/* dot product of row j of Tip and row i of X	*/
			}
			sumsq+=sum*sum;
			*ppV++=sum; /* i,j element of XTi */
		}
		(*pd).i=i;
		(*pd++).d=sumsq;
		if (criterion) { /* row of XTi times Ti' making (XTiTi') */
			pU=U+i*k;
			q=Ti;
			for (j=0;j<k;j++) {
				p=pV++;
				sum=0;
				for (l=j;l<k;l++) {
					sum+=(*p++)*(*q++);
				}
				*(pU++)=sum;
			}
		}
		pX+=k; /* next row of X */
	}
	if (criterion==2) {		
		MatMult(B,U,BU,k,N);
	}
	dShellSort(d,N,true); /* sort d in increasing order */

	*maxd=d[N-1].d;
	*maxdi=d[N-1].i;

	g=0;
	for (i=0;i<N;i++) {
		if (!doApprox && designFlag[d[i].i]) /* designFlag marks nullify pts when doApprox is true */
			d[i].o=g++; /* design points separately ordered */
		else
			d[i].o=i;			/* put order in o */
	}
	dShellSort(d,N,false);	 /* return to original order */
}


/* findDeltaAlpha ***********************************************************************
|	Calculates delta and alpha, for a normalized dispersion matrix
|	For the D crierion it is delta=d, where d=max(d(u)), and alpha=(d-k)/((d-1)k)
|   For a linear criterion L(A) it is delta=phi where phi=max(phi(u)), 
|		and alpha=(phi-L(Mi))/(phi(d(u)-1))
|	For A, phi(u) is trace(u'Mi*Mi u), where Mi is the inverse of the normalized dispersion matrix,
|		and L(Mi)=trace(Mi)	
|	For I, phi(u) is trace(u'Mi*B*Mi,u), and L(Mi)=trace(BMi),
|   See Federov section 2.10 - 2.12.
|	returns delta, xnew and alpha
*/

double findDeltaAlpha(
	double  *bestAlpha,
	double  *BU,		/* Pointer to BU, which is U transformed by whole space crossproduct matrix, G, used by I only */
	int     criterion, /* 0 for D, 1 for A, 2 for I */
	int		*xnew,
	double  maxd,
	int		maxdi,
	double  Acrit,
	double  Icrit,
	dType	*d,
	double	*U, /* U=XM^-1 */
	double	N,
	int		k,
	bool    *failure
)
{
	int		i;
	double	*pUi=0;	/* ith row of U */
	double  *pBUi=0;   /* ith row of BU */
	double	bestDelta=1e-14;
	int		xi=-1;		/* row number for extremeDelta */
	double  crit=(criterion==0)?k:((criterion==1)?Acrit:Icrit);
	double  delta;
	double  alpha;

	*failure=false;


	if (criterion==0) {
		bestDelta=maxd; /* Federov 2.5.11 */
		*bestAlpha=(bestDelta-crit)/((bestDelta-1)*(double)k);
		bestDelta-=crit;
		xi=maxdi;
	}
	else {
		for (i=0;i<N;i++) {
			if (criterion) /* A and I */
				pUi=U+k*i;
			if (criterion==2) /* I */
				pBUi=BU+k*i;

			delta=GetLinearCriterionA(pBUi,criterion,pUi,k);
			alpha=(delta-crit)/(1.2*delta*(d[i].d-1)); /* Federov 2.10.10 */
			delta-=crit;

			if (delta>bestDelta) {
				*bestAlpha=alpha;
				bestDelta=delta;
				xi=i;
			}
		}
	}

	if (xi==-1) {
		*failure=true;  /* fails to find a delta */
	}

	*xnew=xi;
	return(bestDelta);
	
}
/* findDelta ***********************************************************************
|	Calculates delta for a normalized dispersion matrix
|	For the D crierion it is delta=(n*dj-(di*dj-dij^2)-n*di)/n^2, where i is in design and j is not
|   For a linear criterion L(A) it is 
|	{(n-di)L(xnew.xnew')+dij[L(xold.xnew')+L(xnew.xold')]-(n+dj)L(xold.xold')}/(n*(1+delta))
|	For A, L(u,v) is v'Mi*Mi u, where Mi is the inverse of the normalized dispersion matrix.
|	For I, L(u,v) is v'Mi*B*Mi,u.
|   See Federov section 3.3, and my notes
|	returns delta and xold and xnew
*/

double findDelta(
	double  *BU,		/* Pointer to BU, which is U transformed by whole space crossproduct matrix, G, used by I only */
	int     criterion, /* 0 for D, 1 for A, 2 for I */
	int		*xold,
	int		*xnew,
	dType	*d,
	double	*U,
	double  *V,
	double	N,
	int		n,
	int		k,
	bool	*designFlag,
	int		*rows,
	bool    *failure
)
{
	double	di;	/* dx for row in design */
	double	dj;	/* dx for row not in design */
	double	dij; /* dxy for i and j */
	double  dn=(double)n;
	int		i;
	int		j;
	int		l;
	double	*pUi=0;	/* ith row of U */
	double	*pUj=0;	/* jth row of U */
	double  *pVi=0;
	double  *pVj=0;
	double  *pBUi=0;   /* ith row of BU */
	double  *pBUj=0;   /* jth row of BU */
	double	*p;
	double	*q;
	double	delta;
	double  deltaT;
	double	extremeDelta=1e-14;
	int		xi=-1;		/* rows[i] for extremeDelta */
	int		xj=-1;		/* j for extremeDelta */

	*failure=false;

	for (i=0;i<n;i++) {
		if (designFlag[rows[i]]==2 || d[rows[i]].o>Klimit)
			continue; /* omit high d design values -- these are really for the D criterion, but calculating */
						/* and using the proper values for other criteria would gain little and be costly. */
						/* Also omit augmentee points */
		di=d[rows[i]].d;
		pVi=V+k*rows[i];
		if (criterion)
			pUi=U+k*rows[i];
		if (criterion==2)
			pBUi=BU+k*rows[i];
		for (j=0;j<N;j++) {
			if (designFlag[j] || d[j].o<Llimit) /* omit design and low value d's */
				continue;
			dj=d[j].d;
			pVj=V+k*j;
			dij=0;
			p=pVi;
			q=pVj;
			for (l=0;l<k;l++) {
				dij+=(*p++)*(*q++);
			}
			deltaT=delta=(dn*dj-(di*dj-dij*dij)-dn*di)/(dn*dn); 
			if (criterion!=0) {
				pUj=U+k*j;
				if (criterion==2)
					pBUj=BU+k*j;
				delta=GetLinearCriterion(pBUi,pBUj,criterion,pUi,pUj,di,dj,dij,k,dn)/(1+deltaT);
			}
			if (delta>extremeDelta) {
				extremeDelta=delta;
				xi=rows[i];
				xj=j;
			}
		}
	}

	if (xi==-1 || xj==-1) 
		*failure=true;  /* fails to find a delta */

	*xold=xi;
	*xnew=xj;
	return(extremeDelta);
	
}

/* GetLinearCriterionA ****************************************************************
|	Calculates the value for a linear criterion for point i
*/

double GetLinearCriterionA(
	double *pBU,
	int criterion, /* 1 for A 2 for I */
	double *pU,
	int    k
)
{
	double	s=0;
	int		i;

	if (criterion==1) {  /* A criterion */
		for (i=0;i<k;i++) {
			s+=pU[i]*pU[i];
		}
	}
	else {  /* I criterion */
		for (i=0;i<k;i++) { 
			s+=pBU[i]*pU[i];
		}
	}
	return(s);

}


/* GetLinearCriterion ****************************************************************
|	Calculates the value for a linear criterion
|	i referes to a point in the design, j to an outside point
*/

double GetLinearCriterion(
	double *pBUi,
	double *pBUj,
	int criterion, /* 1 for A 2 for I */
	double *pUi,
	double *pUj,
	double di,
	double dj,
	double dij,
	int    k,
	double dn
)
{
	double	si=0;
	double  sji=0;
	double  sij=0;
	double  sj=0;
	double  t=0;
	double  ss;
	int		i;
	double  dn2=dn*dn;

	if (criterion==1) {  /* A criterion */
		for (i=0;i<k;i++) {
			si+=(*pUi)*(*pUi);
			sji+=(*pUi++)*(*pUj);
			t=*pUj++;
			sj+=t*t;
		}
		ss=2*sji;
	}
	else {  /* I criterion */
		for (i=0;i<k;i++) { 
			si+=(*pBUi)*(*pUi);
			sji+=(*pBUj)*(*pUi++);
			sij+=(*pBUi++)*(*pUj);
			sj+=(*pBUj++)*(*pUj++);
		}
		ss=sji+sij;
	}
	return((sj*(dn-di)+ss*dij-si*(dn+dj))/dn2);

}



/* updateA **************************************************************************
|	Updates T and proportions by adding xnew with proportion alpha.
*/
void updateA(
	int		xnew,
	double  *proportions,
	double  alpha,
	double	*T,
	double	*X,
	double	*tVec,
	int		k,
	int		N
)
{
	double	*pXn=X+k*xnew; /* pointer to new row of X */
	double  *pT;
	double  beta=1-alpha;
	int		i;


	pT=T; /* rescale the old design by 1-alpha */
	for (i=0;i<k;i++) {
		(*pT)*=beta;
		pT+=k-i;
	}
		/* Adds new point with weight alpha -- may duplicate an existing point */
	Rotate(pXn,tVec,T,k,k,alpha,1.0);

		/* Adjust the proportions */
	for (i=0;i<N;i++) {
		proportions[i]*=beta;
	}
	proportions[xnew]+=alpha;
}

/* update **************************************************************************
|	Updates T by replacing row xold with row xnew.
|	resets rows and designFlag to conform
*/
void update(
	int		xold,
	int		xnew,
	int		*rows,
	int		*designFlag,
	double	*T,
	double	*X,
	double	*tVec,
	int		k,
	int		n
)
{
	double	*pXo=X+k*xold; /* pointer to old row of X */
	double	*pXn=X+k*xnew; /* pointer to new row of X */
	int		i;


	Rotate(pXn,tVec,T,k,k,1.0,(double)n);
	Rotate(pXo,tVec,T,k,k,-1.0,(double)n);
	designFlag[xold]=false;
	designFlag[xnew]=true;


	for (i=0;i<n;i++) {
		if (xold==rows[i]) {
			rows[i]=xnew;
			break;
		}
	}



	
}




/* FillInB ****************************************************************************
|  Fills in B=upper triangle of X'X/N
*/

void FillInB(
	double *X,
	double *B,
	int k,
	int N
)
{
	double *pB;
	double *pX1;
	double *pX2;
	int i;
	int j;
	int l;
	int K=k*(k+1)/2;
	
	memset((void *)B,0,sizeof(double)*K);
	for (l=0;l<N;l++) {
		pB=B;	   /* Start of B */
		pX1=X+k*l; /* next row of X */
		for (i=0;i<k;i++) {
			pX2=pX1;
			for (j=i;j<k;j++) {
				(*pB++)+=(*pX1)*(*pX2++)/(double)N;		
			}
			pX1++;
		}
	}	
}





/* ProgAlloc ***********************************************************************
|	Allocates space for FederovOptimize() using R_alloc()
|	Change as needed for particular OS
|
|
|	doubles:
|		U:	Nxk
|		V:  Nxk
|		B:  k*(k+1)/2
|		BU: Nxk,
|		T:	k*(k+1)/2
|		Ti: k*(k+1)/2
|		Tip:	k*(k+1)/2
|		W:	k*(k+1)/2
|		vec:	k
|	structure dType
|		d:  N
|	ints
|		designFlag: N
|		trows: m+1
|       
*/


int ProgAlloc(
	double	**U,	
	double  **V, 
	double  **B,   
	double  **BU,   
	double	**T,
	double  **Ti,	
	double	**Tip, 
	double  **W,
	double  **maxmin,
	dType	**d,   
	double  **vec, 
	int		**designFlag, 
	int     **ttrows,
	int		**trows,	
	int		N,  
	int		n, 
	int		k, 
	bool    criterion,
	bool	evaluateI,
	bool    doSpace
)
{
	int K=k*(k+1)/2;

	*V=(double *)R_alloc(N*k,sizeof(double));
	if (!*V) return 3;
	if (criterion) {
		*U=(double *)R_alloc(N*k,sizeof(double));
		if (!*U) return 4;
	}
	if (criterion==2 || evaluateI) {
		if (!doSpace) {
			*B=(double *)R_alloc(K,sizeof(double));
			if (!*B) return 4;
		}
		*BU=(double *)R_alloc(N*k,sizeof(double));
		if (!*BU) return 4;
		*W=(double *)R_alloc(K,sizeof(double));
		if (!*W) return 7;
	}
	*T=(double *)R_alloc(K,sizeof(double));
	if (!*T) return 5;
	*Ti=(double *)R_alloc(K,sizeof(double));
	if (!*Ti) return 5;
	*Tip=(double *)R_alloc(K,sizeof(double));
	if (!*Tip) return 6;
	*maxmin=(double *)R_alloc(2*k,sizeof(double));
	if (!maxmin) return 7;
	*d=(dType *)R_alloc(N,sizeof(dType));
	if (!*d) return 8;
	*vec=(double *)R_alloc(k,sizeof(double));
	if (!*vec) return 9;
	*designFlag=(int *)R_alloc(N,sizeof(int));
	if (!*designFlag) return 10;
	*ttrows=(int *)R_alloc(N,sizeof(int));
	if (!*ttrows) return 10;
	*trows=(int *)R_alloc(N,sizeof(int)); /* N because Permute needs N for RandomStart */
	if (!*trows) return 11;

	return 0;
}

/* ProgDealloc *********************************************************************
|	Frees space allocated by ProgAlloc
*/

void ProgDealloc(
	double	*U,	
	double  *V,  
	double  *B,
	double  *BU,
	double	*T,	
	double  *Ti,  
	double	*Tip,  
	double  *W,	 
	double  *maxmin,
	dType	*d,  
	double  *vec, 
	int 	*designFlag, 
	int		*ttrows,
	int		*trows
)	
{
	if (U)
		Free(U);
	if (V)
		Free(V);
	if (B)
		Free(B);
	if (BU)
		Free(BU);
	if (T)
		Free(T);
	if (Ti)
		Free(Ti);
	if (Tip)
		Free(Tip);
	if (W)
		Free(W);
	if (maxmin)
		Free(maxmin);
	if (d)
		Free(d);
	if (vec)
		Free(vec);
	if (designFlag)
		Free(designFlag);
	if (ttrows)
		Free(ttrows);
	if (trows)
		Free(trows);

}

/* filloutDesign ********************************************************************
|  Adds rows to a nullification start until n is reached
*/
void filloutDesign(
	double	*T,
	double	*Ti,
	double	*Tip,
	int		*rows,
	int		*ttrows,
	double	*X,
	double	*maxmin,
	double	*vec,
	int		k,
	int		ka,
	int		n,
	int     N,
	bool	*singular
  )
{
	int i;
	int j;
	int l;
	int m;
	int newRow;
	double *pTip;
	double *pX;
	double *ppX;
	double delta;
	double sumsq;
	double d;
	double norm=1.0;

	

		/* Note that in this function, T is not scaled by n; hence it must be scaled to */
	    /* be usable elsewhere. At the present time, the return value of T is not used elsewhere. */


	reduceDesign(rows,X,T,maxmin,vec,k,ka,false,singular);

	
	if (*singular) {
		return;
	}

	makeTiAndTipFromT(Tip,T,Ti,maxmin,norm,singular,k);
	if (*singular)
		return;

	for (i=ka;i<n;i++) {  
		newRow=-1;
		delta=-1;
		for (j=0;j<N;j++) {
			sumsq=0;			/* d=(xTiTi'x'), where x is a row of X and then find delta=max(d) */
			if (!ttrows[j]) {
				pX=X+j*k;
				pTip=Tip;
				for (l=0;l<k;l++) {
					ppX=pX;
					d=0;
					m=l+1;
					while (m--) {
						d+=(*ppX++)*(*pTip++);
					}
					sumsq+=d*d;
				}
				if (sumsq>delta) {
					delta=sumsq;
					newRow=j;
				}
			}			

		}
		if (newRow==-1) {
			*singular=true;
			return; /* singular */
		} else {
			ttrows[newRow]=1;
		}

		pX=X+newRow*k;
		rows[i]=newRow;

		if (i!=n-1) {
				/* Add the row to the design */
			Rotate(pX,vec,T,k,k,1.0,1.0);  
			makeTiAndTipFromT(Tip,T,Ti,maxmin,norm,singular,k);
			if (*singular)
				return;
		}
	}

	*singular=false;

}

/* FederovOptimize ********************************************************************
|	The optimizing function. 
|	Returns 12 in case the determinant is not positive
|	otherwise returns the number of iterations.
|
*/

#define designTol 1e-3

int FederovOptimize(
	double	*X, /* candidate list of dim Nxk, where k is the number of model terms */
	double  *B, /* Upper triangle of Space matrix, G=X'X/N, for I or NULL if criterion!=2 */
	double  *BU, /* Product of G and U or NULL if criterion !=2 */
	double  *proportions, /* Proportions for approximate theory or NULL if not done */
	bool	RandomStart, /* if TRUE rows will be filled at random to start */
	int 	Nullify,   /* if non-zwro nullify() will be used to fill in rows: if 2, some pts will be random */
	int		criterion,  /* 0 for D, 1 for A, 2 for I */
	bool	evaluateI, /* Return I as well as D and A */
	bool    doSpace,
	int		augment,   /* number of rows for aumentee design*/
	double  *D, /* D determinant always returned*/
	double	*A, /* A  always returned */
	double	*I, /* I eff returned only when criterion==2 or evaluateI is true */
	double  *G, /* G eff always returned*/
	double	*U,	/* U=XTiTi', Nxk or NULL if criterion!=0 */
	double  *V, /* V=XTi Nxk */
	double	*T,	/* Design reduced, to T,upper triangular kxk */
	double  *Ti, /* T^-1 */
	double	*Tip,  /* Ti' rearranged so column i of Ti becomes row i of Tip */
	double  *W,	 /* working storage */
	double  *maxmin, /* maximum and minumum ranges. */
	dType	*d,   /* (d_i=U_i*U_i) */
	double  *vec, /* Temporary storage */
	int 	*designFlag, /* non-zero if row of X is in design Set to 2 for augmentee points */
						 /* It is used to mark nullfy induced non support points when doAapprox is true. */
	int		*ttrows, /* Used only in filloutDesign() to avoid duplicate rows when nullify=1 */
	int		*rows,	/* input/output array of indexes of design row numbers from X */
	int		*trows,	/* working array */
	int     N,  /* Number of rows in X */
	int		n,	/* Number of rows in design */
	int		k,	    /* Number of model terms */
	int		maxIteration,
	int     nRepeats, 
	double  DFrac,	/* fraction of design points included in search */
	double  CFrac, /* fraction of candidate points included in search */
	int     *error
)
{

	double	delta=0;
	double	maxd;
	int		maxdi;
	int		xold,
			xnew;
	int		iter=0;
	int		miter;
	double  crit;
//	double  logDcrit;
	double  Acrit;
	double  Icrit=0;
	double  logdet;
	bool	singular;
	int		countSingular=0;
	int		nRepeatCounts=nRepeats;
	double  bestCrit=criterion?1e8:-1e8;
	int		i;
	int     j;
	int     l;
	int     pp;
	int		ka;
	double  norm=(doApprox)?1:(double)n;
	bool	failure;
	double  alpha;
	double  pf;
	bool	firstPass=true;
	int		nd;
	int		np;



	if (!doSpace && (criterion==2 || evaluateI)) 
		FillInB(X,B,k,N);

					/* Recommend KFrac=1, LFrac=0 */
					/* Use KFrac=0, LFrac=1 to use the two extreme values only */
	Klimit=(int)(DFrac*n); /* uses only d values less than Klimit -- if zero, only the smallest is used. */
	Llimit=(int)((1-CFrac)*(N-1)); /* uses only d values greater than Llimit -- if N-1, only the largest is used. */


	if (maxIteration==0)  
		maxIteration=100;

	if (doApprox) {
		n=k;
		Nullify=1; /* Never start with a random design, since it will too often by */
				   /* singular when n is exactly k */
	}



	repeat{

		firstPass=true;

		memset((void *)designFlag,0,N*sizeof(int));
		if (augment) {
			for (i=0;i<augment;i++) {
				designFlag[rows[i]]=2; /* These rows will never be exchanged */
			}
		}

		if (Nullify) {

			if (augment<k) {  /* Finds k starting rows */
				if (!nullify(X,V,augment,rows,designFlag,N,k)) {
					*error=13;
					return(0);
				}
			}
				
			for (i=0;i<N;i++) {
				ttrows[i]=designFlag[i];
			}

			ka=maxm(k,augment);

			if (n>ka) {	/* Adds rows to make n by choosing points with max variance */
				if (Nullify==1) {
					filloutDesign(T,Ti,Tip,rows,ttrows,X,maxmin,vec,k,ka,n,N,&singular);
					if (singular) {
						countSingular++;
						break;
					}
					for (i=0;i<n;i++)
						trows[i]=rows[i];
				}
				else {	/* Adds rows to make n at random */
					j=0;
					for (i=0;i<N;i++) {
						if (!designFlag[i])
							ttrows[j++]=i;
					}
					/* add points at random to fill out the design.  */
					Permute(ttrows,N-k);
					j=0;
					l=0;
					for (i=0;i<N;i++) {
						if (designFlag[i]) {
							pp=rows[j];
							trows[j++]=pp;
						}
						else {
							pp=ttrows[l];
							trows[l++]=pp;
						}
					}
/*
					j=n-k;
					for (i=0;i<k;i++) {
						trows[j++]=rows[i];
					}
*/
				}
			}
			else {
				for(i=0;i<n;i++)
					trows[i]=rows[i];
			}
		} else
		if (RandomStart) {			/* Just a random selection of rows */
			for (i=0;i<N;i++) {
				trows[i]=i;
			}
			Permute(trows,N);
		} else {
			if (augment!=0) {		/* Augmenting. Add rows at random to the augment rows */
				for (i=0;i<augment;i++) {
					trows[i]=rows[i];
				}
				if (augment<n) {
					j=augment;
					for (i=0;i<N;i++) {
						if (!designFlag[i]) {
							trows[j++]=i;
						}
					}
					Permute(trows+augment,N-augment);
				}
			}
			else {					/* The user has supplied the rows */
				for (i=0;i<n;i++) {
					trows[i]=rows[i];
				}
			}
		}

        for (i=0;i<n;i++) {
			if (!designFlag[trows[i]])
				designFlag[trows[i]]=1; /* designFlag will have 2 for augment rows and 1 for others that may */
										/* be exchanged. */
		}

		/* The T from nullify() is not scaled by n; hence remake it here. */
        reduceDesign(trows,X,T,maxmin,vec,k,n,true,&singular);

		if (doApprox) { /* Set the initial proportions to 1/k for rows in the design. */
			memset((void *)proportions,0,N*sizeof(double));
			pf=1.0/(double)k;
			for (i=0;i<k;i++) {
				proportions[rows[i]]=pf;
			}
		}

restart: if (!singular) {
			logdet=makeTiAndTipFromT(Tip,T,Ti,maxmin,norm,&singular,k);
			if (singular)
				goto sing;
			crit=evaluateCriteria(Tip,Ti,W,B,criterion,evaluateI,&Acrit,&Icrit,logdet,k,N);
			makeBXd(X,U,V,Ti,Tip,B,BU,criterion,designFlag,d,&maxd,&maxdi,k,N);
			miter=0;
			repeat
				if (doApprox) {
					delta=findDeltaAlpha(&alpha,BU,criterion,&xnew,maxd,maxdi,Acrit,Icrit,d,U,N,k,&failure);
					if (failure || delta<((crit>0)?crit*designTol:designTol) || miter>maxIteration) {
						break;
					}
					R_CheckUserInterrupt();
					updateA(xnew,proportions,alpha,T,X,vec,k,N);
					designFlag[xnew]=0; /* The point is a support point, so remove it from the nullify */
										/* induced points */

				}
				else {
					delta=findDelta(BU,criterion,&xold,&xnew,d,U,V,N,n,k,designFlag,trows,&failure);
					if (failure || delta<crit*designTol|| miter>maxIteration) {
						break;
					}
					R_CheckUserInterrupt();
					update(xold,xnew,trows,designFlag,T,X,vec,k,n);
				}
				logdet=makeTiAndTipFromT(Tip,T,Ti,maxmin,norm,&singular,k);
				if (singular)
					goto sing;
				crit=evaluateCriteria(Tip,Ti,W,B,criterion,evaluateI,&Acrit,&Icrit,logdet,k,N);
				makeBXd(X,U,V,Ti,Tip,B,BU,criterion,designFlag,d,&maxd,&maxdi,k,N);
				miter+=1;
			forever;
			if (miter>iter)
				iter=miter;
		
			if (criterion?crit<bestCrit:crit>bestCrit) {
				bestCrit=crit;
				if (!doApprox) {
					for (i=0;i<n;i++)
						rows[i]=trows[i];
				}
				*D=exp(logdet/(double)k);
				*A=Acrit;
				*G=k/(maxd);
				if (criterion==2 || evaluateI)
					*I=Icrit;
			}

			if (doApprox && firstPass) { /* Remove nullify induced points and run again */
				firstPass=false;
				np=0;
				nd=0;
				for (i=0;i<N;i++) {
					if (proportions[i]>0)
						np++;
					if (designFlag[i]>0)
						nd++;
				}
				n=np-nd;
				if (n>=k) { /* Are there enough non nullify points to support the model? */
					for (i=0;i<N;i++) {
						if (designFlag[i]>0)
							proportions[i]=0;
					}
					memset((void *)designFlag,0,N*sizeof(int)); /* not used again */
					j=0;
					for (i=0;i<N;i++) {
						if (proportions[i]>0) {
							trows[j++]=i;
							proportions[i]=1/(double)n;
						}
					}
					reduceDesign(trows,X,T,maxmin,vec,k,n,true,&singular);
					goto restart;
				}
			}
		}
		else {
sing:			countSingular++;
		}
	}until((!RandomStart && Nullify!=2) || !(--nRepeatCounts)); 

	if (RandomStart || Nullify==2) {
		if (countSingular==nRepeats)
			*error=12;
	}
	else {
		if (countSingular)
			*error=12;
	}


	return(iter);

}


/* transposeMatrix ****************************************************************************
|  Transposes a matrix in place
|  I got this from Robin Becker robin@jessikat.fsnet.co.uk
*/

void transposeMatrix(
	double *X,
	int N,  /* number of rows in transposed matrix */
	int k	/* number of columns in transposed matrix */
)	
{
    int     size = N*k-2;
	int		i=1;
	int		row;
	int		column;
	int		indx;
	double  temp;

	for(i=1;i<size;i++){
		indx = i;
		repeat
			column = indx/k; /* indx (column+k*row) of cell in transposed matrix */
			row = indx%k;
			indx = N*row +  column; /* the value to be swapped is in cell at indx in current matrix */
		solongas(indx < i); /* indx<i is a previously exchanged cell, */
							/* track its original contents and then exchange when found */

		if (indx >i) {
			temp=X[i];
			X[i] = X[indx];
			X[indx] = temp;
		}
	}
}



/* FederovOpt *****************************************************************************
|	calling function for FederovOptimize() with memory allocation and deallocation
|
|
|	Input: 
|		X, an matrix of model expanded candidate points
|		RandomStart, if TRUE rows will be filled at random to start. If FALSE, the input
|			rows will be used to start.
|		Nullify. If non-zero, nullification will find starting design -- if 2 some points will be random 
|		rows, a vector of row numbers from X, representing the starting design -- space
|			must always be allocated in the calling program. If doApprox is TRUE, rows 
|			must be at least k.
|		criterion, the criterion (0 for D, 1 for A, 2 for I)
|		evaluateI, if true I will be evaluated, otherwise only D and A
|		doSpace,  when true the user has input B
|		B, the space matrix. upper triangular part of S'S/ncol(S)
|		augment, the number of rows in Xin that should not be exchanged
|		proportions, a vector of N doubles or NULL if approximate theory is not done
|		n, number of rows in the design -- ignored if doApprox is TRUE
|		maxIteration, maximum iterations within algorithm
|		nRepeats, number of retries -- each from a different random start 
|		DFrac, fraction of smallest variance design points included in search
|		CFrac, fraction of smallest variance candidate points included insearch
|	Output: 
|			function return: the number of iterations used
|			rows will contain row numbers of X of the design if doApprox is FALSE
|			proportions will contain the proportions for an approximate theory design
|			D criterion value = (det(M))^(1/k), where M is normalized cross product matrix
|			A criterion value = (trace(MI))/k, where MI is the inverse M
|			I criterion value = trace(G*MI), where G=X'X/N 
|			G efficincy		  = k/max(d), where d is vector of normalized variances in space
|			iter max number of iterations
|			error true when no positive determinant can be found
|	The paramaters retries, maxIteration, KFrac, and LFrac control the program's behavior:
|		retries  mxIteration KFrac LFac
|			5(0)	100(0)	   0.3  0.9   if retries and maxIteration are 0
|											program sets them to 5 and 100
|			5(0)	100(0)	   1     0    normal choice, all points included in search
|			1       anything   0     1    Only lowest and highest single point included in
|										    search, and only one pass is run
| ***********************
|	Allocations and deallocations are handled by ProgAlloc() and ProgDealloc()
|
|	User allocations:
|	double
|		X: Nxk, candidate list
|		rows: n, number of rows in the design, or k when doing approximate theory.
|			These are row numbers from X. It must be of length N.
|		proportions: a vector of N proportions or NULL if approx theory not done
|
|	Program allocations:
|	doubles:
|		U:  Nxk, if (criterion!=0)
|		V:  Nxk.
|		B: k*(k+1)/2, Null or the upper triangular portion of G=X'X a symmetrix crossproduct matrix over all 
|						points in the space. if (criterion==2)
|		BU: Nxk, Null unless  (criterion==2)
|		T:	k*(k+1)/2
|		Tip:	k*(k+1)/2
|		W:	k*(k+1)/2
|		maxmin: 2k
|		vec:	k
|	structure dType
|		d:	N
|	ints 
|		designFlag: N
|		trows: N
*/



SEXP FederovOpt(
	SEXP Xi,
	SEXP RandomStarti,
	SEXP rowsi,
	SEXP Nullifyi,
	SEXP criterioni,
	SEXP evaluateIi,
	SEXP doSpacei,
	SEXP Bi,
	SEXP augmenti,
	SEXP doApproxi,
	SEXP proportionsi,
	SEXP ni,
	SEXP maxIterationi,
	SEXP nRepeatsi,
	SEXP DFraci,
	SEXP CFraci
)
{
	double	*X;
	int		N;
	int		k;
	int		n;
	bool	RandomStart;
	int		*rows;
	int		Nullify;
	int		criterion;
	bool	evaluateI;
	int		augment;
	double	*proportions=0;
	int		maxIteration;
	int		nRepeats;
	double	DFrac;
    double	CFrac;
	int		iter;
	int		error;
		/* returns */
	SEXP	alist=0;
	SEXP    anames;
	SEXP    DVector;
	SEXP    AVector;
	SEXP    IVector;
	SEXP    GVector;
	SEXP    iterVector;
	SEXP    errVector;
	SEXP    rowVector;
	SEXP    propVector;

	double  D;
	double  A;
	double  I;
	double  G;

	bool    doSpace=false;
	double  *B=NULL;
	double  *BU=0;
	double	*U=0;
	double  *V=0;
	double	*T=0;	
	double  *Ti=0;
	double	*Tip=0; 
	double  *maxmin;
	double  *W=0;	
	dType	*d=0; 
	double  *vec=0;
	int 	*designFlag=0; 
	int		*ttrows=0;
	int		*trows=0;
	int     i;



	PROTECT(Xi=AS_NUMERIC(Xi)); /* The only argument modified by the C program */
	X=NUMERIC_POINTER(Xi);

	rows=INTEGER_POINTER(rowsi);
	doApprox=INTEGER_POINTER(doApproxi)[0]; /* A global variable */
	if (doApprox) {
		proportions=NUMERIC_POINTER(proportionsi);
	}

	doSpace=INTEGER_POINTER(doSpacei)[0];
	if (doSpace) {
	    PROTECT(Bi=AS_NUMERIC(Bi));
	    B=NUMERIC_POINTER(Bi);
	}

	N=INTEGER_POINTER(GET_DIM(Xi))[0];
	k=INTEGER_POINTER(GET_DIM(Xi))[1];
	RandomStart=INTEGER_POINTER(RandomStarti)[0];
	Nullify=INTEGER_POINTER(Nullifyi)[0];
	criterion=INTEGER_POINTER(criterioni)[0];
	evaluateI=INTEGER_POINTER(evaluateIi)[0];
	augment=INTEGER_POINTER(augmenti)[0];
	n=INTEGER_POINTER(ni)[0];
	maxIteration=INTEGER_POINTER(maxIterationi)[0];
	nRepeats=INTEGER_POINTER(nRepeatsi)[0];
	DFrac=NUMERIC_POINTER(DFraci)[0];
	CFrac=NUMERIC_POINTER(CFraci)[0];

	if (doApprox || Nullify==1)
		nRepeats=1; /* No random start since more repeats make no sense. */



	transposeMatrix(X,N,k);

	if (!(error=ProgAlloc(&U,&V,&B,&BU,&T,&Ti,&Tip,&W,&maxmin,&d,&vec,
			&designFlag,&ttrows,&trows,N,n,k,criterion,evaluateI,doSpace))) {

		iter=FederovOptimize(X,B,BU,proportions, 
			RandomStart,Nullify,criterion,evaluateI,doSpace,augment,&D,
			&A,&I,&G,U,V,T,Ti,Tip,W,maxmin,d,vec,designFlag,ttrows,rows,trows,N,n,k,
			maxIteration,nRepeats,DFrac,CFrac,&error);

			/* Prepare a list to return */
		PROTECT(alist=NEW_LIST(8));

		PROTECT(DVector=NEW_NUMERIC(1));
		NUMERIC_POINTER(DVector)[0]=D;
		SET_ELEMENT(alist,0,DVector);
		UNPROTECT(1);

		PROTECT(AVector=NEW_NUMERIC(1));
		NUMERIC_POINTER(AVector)[0]=A;
		SET_ELEMENT(alist,1,AVector);
		UNPROTECT(1);

		PROTECT(IVector=NEW_NUMERIC(1));
		NUMERIC_POINTER(IVector)[0]=I;
		SET_ELEMENT(alist,2,IVector);
		UNPROTECT(1);

		PROTECT(GVector=NEW_NUMERIC(1));
		NUMERIC_POINTER(GVector)[0]=G;
		SET_ELEMENT(alist,3,GVector);
		UNPROTECT(1);

		PROTECT(iterVector=NEW_INTEGER(1));
		INTEGER_POINTER(iterVector)[0]=iter;
		SET_ELEMENT(alist,4,iterVector);
		UNPROTECT(1);

		PROTECT(errVector=NEW_INTEGER(1));
		INTEGER_POINTER(errVector)[0]=error;
		SET_ELEMENT(alist,5,errVector);
		UNPROTECT(1);

		// CRAN correction 2014-02-11
		int Nr = LENGTH(rowsi);
		PROTECT(rowVector=NEW_INTEGER(Nr));
		for (i=0;i<Nr;i++) {
			INTEGER_POINTER(rowVector)[i]=rows[i];
		}
		SET_ELEMENT(alist,6,rowVector);
		UNPROTECT(1);

		if (doApprox) {
			PROTECT(propVector=NEW_NUMERIC(N));
			for (i=0;i<N;i++) {
				NUMERIC_POINTER(propVector)[i]=proportions[i];
			}
			SET_ELEMENT(alist,7,propVector);
			UNPROTECT(1);
		}

			/* Label the variables in the list */

		PROTECT(anames=NEW_CHARACTER(8));
		SET_STRING_ELT(anames,0,mkChar("D"));
		SET_STRING_ELT(anames,1,mkChar("A"));
		SET_STRING_ELT(anames,2,mkChar("I"));
		SET_STRING_ELT(anames,3,mkChar("G"));
		SET_STRING_ELT(anames,4,mkChar("iter"));
		SET_STRING_ELT(anames,5,mkChar("error"));
		SET_STRING_ELT(anames,6,mkChar("rows"));
		SET_STRING_ELT(anames,7,mkChar("proportions"));
		SET_NAMES(alist,anames);
		UNPROTECT(1);

		UNPROTECT(1);

		if (doSpace)
		    UNPROTECT(1); /* Bi */

		UNPROTECT(1);	/* Xi */

		return alist;
	}
	else {
		if (doSpace)
		    UNPROTECT(1); /* Bi */

		UNPROTECT(1);	/* Xi */

		return R_NilValue;
	}

	/* In case Calloc() is used */
/*	ProgDealloc(U,V,B,BU,T,Ti,Tip,W,maxmin,d,vec,designFlag,ttrows,trows); */


}
