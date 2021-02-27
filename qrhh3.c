#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bla.h"
#include "qrhh3.h"

void qrhh3(int m, double **A) {
	/* Compute a QR factorization of the square upper Hessenberg
	 * matrix A, by Householder reflections.  A=QR.   A is 
	 * overwritten by RQ.  This routine takes advantage of the
	 * fact that A is upper Hessenberg, which makes the reflection
	 * vectors v have at most two nonzero entries. */
	int i,j,k;
	double **v,normv,sx1,vw;

	/* Because H is upper Hessenberg, to zero out below the diagonal
	 * in column k, the vector v[k] onto which we need to project will
	 * only have two nonzero entries in general.   We store just these 
	 * two entries as a row in a matrix called v.  Each of these rows
	 * are normalized to length one.  If we wish to indicate that no
	 * reflection needs to be done to this column (because the 
	 * subdiagonal entry is already zero) then we store a 3.0 in the 
	 * first entry of v[k], and we test via   if (v[k][0] < 2.0).  */
	if ((v = (double **) malloc((m-1)*sizeof(double *))) == NULL) {
		fprintf(stderr,"malloc failed.\n");
		exit(1);
	}
	for (k=0; k<m-1; k++) {
		if ((v[k] = (double *) malloc(2*sizeof(double))) == NULL) {
			fprintf(stderr,"malloc failed.\n");
			exit(1);
		}
	}
	for (k=0; k<m-1; k++) {
		/* If the subdiagonal entry is already zero, then nothing to do
		 * this column.  Just store a 3.0 in v[k][0] to signal this. */
		if (A[k+1][k] == 0.0)
			v[k][0] = 3.0;
		else {
			/* Determine the vector v that x will be projected onto. */
			for (i=0; i<2; i++) v[k][i] = A[i+k][k];
			normv = norm(2,v[k]);
			if (v[k][0] >= 0.0)
				sx1 = normv;
			else
				sx1 = -normv;
			v[k][0] += sx1;
			/* Normalize v. */
			normv = normalize(2,v[k]);
			/* Apply the reflection (subtract twice the projection onto v) to 
			 * each column of A from column k to end, rows k to end.   We know
			 * that the kth column projects to [-sx1,0,0,...,0]^T. */
			A[k][k] = -sx1;
			if (k<m-1) A[k+1][k] = 0.0;
			for (j=k+1; j<m; j++) {
				/* Get inner product of column j of A rows k and k+1 with v. */
				vw = inner_col(2,j,A + k,v[k]);
				/* Update rows k and k+1. */
				for (i=0; i<2; i++) A[i+k][j] -= 2*v[k][i]*vw;
			}
		}
	}
	/* A currently stores R.  Overwrite it with R*Q, where Q is defined by
	 * the sequence of vectors v that we have stored. */
	for (k=0; k<m-1; k++) {
		/* Apply the kth reflection on the right to each row of R, only columns
		 * k and k+1 are affected. */
		if (v[k][0] < 2.0) {  // Reflector used for this column. 
			for (i=0; i<m; i++) {
				/* Get inner product of row i of R columns k and k+1 with v. */
				vw = inner(2,A[i]+k,v[k]);
				/* Update columns k and k+1. */
				for (j=0; j<2; j++) A[i][j+k] -= 2*v[k][j]*vw;
			}
		}
	}
	for (k=0; k<m-1; k++) free(v[k]);
	free(v);
	return;
}
