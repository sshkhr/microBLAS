#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "bla.h"
#include "forbacksubs.h"

double *backsubs(int n, double **R, double *b) {
	/* Apply back substitution for the square upper triangular system Rx=b. */
	/* INPUT:
	 * n : number of columns
	 * R : (assumed square) matrix (only upper triangle is accessed)
	 * b : vector
	 *
	 * RETURNS:
	 * x : pointer to solution, NULL if no solution
	 */
	int  i;
	double *x, temp;

	/* Allocate space for the solution vector. */
	if ((x = (double *) malloc(n*sizeof(double))) == NULL) {
		fprintf(stderr,"backsubs: malloc failed\n");
		exit(1);
	}
	/* Work through rows backwards from bottom. */
	for (i=n-1; i>=0; i--) {
		temp = b[i] - inner(n-1-i,R[i]+(i+1),x+(i+1));
		if (R[i][i] == 0.0) {
			/* Warn if no solution, return NULL. */
			if (fabs(temp) > fabs(b[i])*1.0e-15) {
				printf("Warning: System has no solution.\n");
				free(x);
				return NULL;
			}
			else x[i] = 0.0;
		}
		else 
			x[i] = temp/R[i][i];
	}
	return x;
}

double *forsubs(int n, double **L, double *b) {
	/* Apply forward substitution for the square lower triangular system Lx=b. */
	/* INPUT:
	 * n : number of columns
	 * L : (assumed square) matrix  (only lower triangle is accessed)
	 * b : vector
	 *
	 * RETURNS:
	 * x : pointer to solution, NULL if no solution
	 */
	int i;
	double *x, temp;

	/* Allocate space for the solution vector. */
	if ((x = (double *) malloc(n*sizeof(double))) == NULL) {
		fprintf(stderr,"forsubs: malloc failed\n");
		exit(1);
	}
	/* Work through rows forwards from top. */
	for (i=0; i<n; i++) {
		temp = b[i] - inner(i,L[i],x);
		if (L[i][i] == 0.0) {
			/* Warn if no solution, return NULL. */
			if (fabs(temp) > fabs(b[i])*1.0e-15) {
				printf("Warning: System has no solution.\n");
				free(x);
				return NULL;
			}
			else x[i] = 0.0;
		}
		else 
			x[i] = temp/L[i][i];
	}
	return x;
}

double *backsubsT(int n, double **RT, double *b) {
	/* Apply back substitution for the square upper triangular system Rx=b,
	 * where R is stored as its transpose, RT. */
	/* INPUT:
	 * n : number of columns
	 * RT : transpose of matrix (only lower triangle of RT is accessed)
	 * b : vector
	 *
	 * RETURNS:
	 * x : pointer to solution, NULL if no solution
	 */
	int  i;
	double *x, temp;

	/* Allocate space for the solution vector. */
	if ((x = (double *) malloc(n*sizeof(double))) == NULL) {
		fprintf(stderr,"backsubsT: malloc failed\n");
		exit(1);
	}
	/* Work through rows backwards from bottom. */
	for (i=n-1; i>=0; i--) {
		temp = b[i] - inner_col(n-1-i,i,RT + i+1),x+(i+1));
		if (RT[i][i] == 0.0) {
			/* Warn if no solution, return NULL. */
			if (fabs(temp) > fabs(b[i])*1.0e-15) {
				printf("Warning: System has no solution.\n");
				free(x);
				return NULL;
			}
			else x[i] = 0.0;
		}
		else 
			x[i] = temp/RT[i][i];
	}
	return x;
}

double *forsubsT(int n, double **LT, double *b) {
	/* Apply forward substitution for the lower triangular system Lx=b,
	 * where L is stored as its transpose, LT. */
	/* INPUT:
	 * n : number of columns
	 * LT : transpose of matrix  (only upper triangle of LT is accessed)
	 * b : vector
	 *
	 * RETURNS:
	 * x : pointer to solution, NULL if no solution
	 */
	int i;
	double *x, temp;

	/* Allocate space for the solution vector. */
	if ((x = (double *) malloc(n*sizeof(double))) == NULL) {
		fprintf(stderr,"forsubsT: malloc failed\n");
		exit(1);
	}
	/* Work through rows forwards from top. */
	for (i=0; i<n; i++) {
		temp = b[i] - inner_col(i,i,LT,x);
		if (LT[i][i] == 0.0) {
			/* Warn if no solution, return NULL. */
			if (fabs(temp) > fabs(b[i])*1.0e-15) {
				printf("Warning: System has no solution.\n");
				free(x);
				return NULL;
			}
			else x[i] = 0.0;
		}
		else 
			x[i] = temp/LT[i][i];
	}
	return x;
}