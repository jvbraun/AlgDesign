/* file Utility.c
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



SEXP	GetFactorial(SEXP Xi,SEXP levelsi,SEXP centri,SEXP factorVeci);

#define Idx(A,B) ((B)+(A)*N)



/* GetFactorial ************************************************************************
| Generates a factorial candidate set according to levels. Keeps matrix in R's row major form.
*/


SEXP GetFactorial(
	SEXP Xi,
	SEXP levelsi,
	SEXP centeri,
	SEXP factorVeci /* a vector with 1's representing factor positions */
)
{
	double *X;
	int *levels;
	int N;
	int nVars;
	bool center;
	int *factorVec;
	int i;
	int j;
	int l;
	double *pX;
	double mid;
	double v;

	Xi=AS_NUMERIC(Xi);
	X=NUMERIC_POINTER(Xi);
	levels=INTEGER_POINTER(levelsi);
	N=INTEGER_POINTER(GET_DIM(Xi))[0];
	nVars=INTEGER_POINTER(GET_DIM(Xi))[1];
	center=INTEGER_POINTER(centeri)[0];
	factorVec=INTEGER_POINTER(factorVeci);


	for (i=0;i<N;i++) {
		j=i;
		l=0;
		for (l=0;l<nVars;l++) {
		 X[Idx(l,i)]=1+j%levels[l];
		 j/=levels[l];
		}
	}

	if (center) {
		pX=X;
		for (i=0;i<nVars;i++) {
			if (!factorVec[i]) {
				mid=(1+levels[i])/2.0;
				if (levels[i]%2) {
					for (j=0;j<N;j++) {
						v=(*pX-mid);
						(*pX++)=v;
					}
				}
				else {
					for (j=0;j<N;j++) {
						v=(*pX-mid)*2.0;
						(*pX++)=v;
					}
				}
			}
			else {
				pX+=N;
			}
			
		}
	}

	return R_NilValue;
}

/* NextCombination *****************************************************************
| This is a version of NEXCOM on p50 of Nijenhuis, A. and Wilf, H.S. (1978).
|	Combinatorial Algorithms. Academic Press.
|	Returns true if this is not the last combination.
*/

void NextCombination(
	int *R,
	int n,
	int k,
	bool *more
)
{
	static int t;
	static int h;

	if (!*more) { /* initialize */
		memset((void *)R,0,k*sizeof(int));
		R[0]=n;
		t=n;
		h=-1;
		*more=(R[k-1]!=n);
		return;		
	}
	if (t>1) h=-1;
	h++;
	t=R[h];
	R[h]=0;;
	R[0]=t-1;
	R[h+1]++;
	*more=(R[k-1]!=n);

}

/* GetMixture ***********************************************************************
| Generates a mixture candidate set. Keeps matrix in R's row major form
*/

SEXP GetMixture(
	SEXP Xi,
	SEXP degreei
)
{
	double *X;
	int N;
	int k;
	int	degree;
	int *R;
	bool more=false;
	int ind;
	int i;

	Xi=AS_NUMERIC(Xi);
	X=NUMERIC_POINTER(Xi);
	degree=INTEGER_POINTER(degreei)[0];
	N=INTEGER_POINTER(GET_DIM(Xi))[0];
	k=INTEGER_POINTER(GET_DIM(Xi))[1];

	R=(int *)R_alloc(k,sizeof(int));


	ind=0;
	repeat {
		NextCombination(R,degree,k,&more);
		for (i=0;i<k;i++) {
			X[Idx(i,ind)]=R[i]/(double)degree;
		}
		ind++;
	}until(!more);

	return R_NilValue;

}

