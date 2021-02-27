#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bla.h"
#include "matvec_read.h"

#include "helper.c"

#define MIN(a, b) ((a) <= (b) ? (a) : (b))

double** bidiag(int m, int n, double** A) 
{
	/* Compute an upper Hessenberg factorization of the square matrix A, by
	* similarity transforms using Householder reflections.
	* A = Q H Q^{-1}
	* H is upper Hessenberg and overwrites A.
	* Q is not returned. */
	int i,j,k;
	double u[m],normu,sx1,uw;
	double v[m],normv,sx2,vw;

	int dim = MIN(n,m);

	double** A_bidiag; 

	A_bidiag = calloc(2,sizeof(double *));
	for (i=0; i<2; i++)
	    A_bidiag[i] = calloc(dim,sizeof(double));

	for (k=0; k<n; k++) // iterate over columns
	{
		// PERFORM LEFT REFLECTION

		for (i=0; i<m-k; i++) 
			u[i] = A[i+k][k];

		// If the second entry down of u are all zeros, nothing to do this column.
		normu = inner(m-k-1,u+1,u+1);
		
		if (normu != 0.0) 
		{
			normu = sqrt(normu + u[0]*u[0]);
			
			if (u[0] >= 0.0)
				sx1 = normu;
			else
				sx1 = -normu;

			u[0] += sx1;
			// Normalize u.
			normu = normalize(m-k,u);

			//printf("Left vector U: ");
			//print_vector(m-k,u);
			
			// Apply the reflection (subtract twice the projection onto v) to
			// each column of A, rows 0 to k are unaffected. Columns 0 to k-1 are
			// unaffected since they have zeros in rows k+1 onward. Rows k+1 to end
			// of column column k project to [-sx1,0,0,...,0]^T.
			
			A[k][k] = -sx1;

			A_bidiag[1][k] = A[k][k]; // store diagonal element in A_bidiag
			
			for (i=k+1; i<m; i++) A[i][k] = 0.0;
			
			for (j=k+1; j<m; j++) 
			{
				// Get inner product of column j of A from row k+1 down with v.
				uw = inner_col(m-k-1,j,A+k+1,v);
				// Update all rows from k+1 down.
				for (i=0; i<m-k-1; i++) 
					A[i+k+1][j] -= 2*u[i]*uw;
			}

			//printf("After Left reflection #%d \n", k);
			//print_matrix(5, 5, A);

		}

		// PERFORM RIGHT REFLECTION
		for (i=0; i<n-k; i++) 
			v[i] = A[k][i+k+1];

		// If the second entry down of u are all zeros, nothing to do this column.
		normv = inner(n-k-2,v+1,v+1);

		if (normv != 0.0) 
		{
			normv = sqrt(normu + v[0]*v[0]);
			
			if (v[0] >= 0.0)
				sx2 = normv;
			else
				sx2 = -normv;

			v[0] += sx2;

			// Normalize v.
			normv = normalize(n-k-1,v);
			
			//printf("Right vector V: ");
			//print_vector(n-k,u);

			// Apply the reflection (subtract twice the projection onto v) to
			// each row of A, rows 0 to k are unaffected. Columns 0 to k-1 are
			// unaffected since they have zeros in rows k+1 onward. Rows k+1 to end
			// of column column k project to [-sx1,0,0,...,0]^T.
			
			A[k][k+1] = -sx2;

			A_bidiag[0][k] = A[k][k+1]; // store super diagonal element in A_bidiag

			for (i=k+2; i<m; i++) A[k][i] = 0.0;

			// Now apply the reflection on the right to each row of A below row k, columns k to n
			// are unaffected.
			for (i=k+1; i<m; i++) 
			{
				// Get inner product of row i of A from column k+1 right with v.
				vw = inner(n-k-1,A[i]+k+1,v);
				
				// Update all columns from k+1 right.
				for (j=0; j<n-k-1; j++) 
					A[i][j+k+1] -= 2*v[j]*vw;
			}

			//printf("After Right reflection #%d \n", k);
			//print_matrix(5, 5, A);

			
		}

		// Store u in column k from row k onwards
		// Store v in row k from column k+1 onwards

		for (i=0; i<m-k; i++) 
			A[i+k][k] = u[i];

		for (i=0; i<n-k; i++) 
			 A[k][i+k+1] = v[i];

	}
	A_bidiag[0][m-2] = A[m-2][m-1]; // store super diagonal element in A_bidiag
	A_bidiag[1][m-1] = A[m-1][m-1]; // store super diagonal element in A_bidiag


	// if number of rows greater than columns fill them with zero
	for (k=n; k<m; k++) // iterate over rows
		for (i=0; i<n; i++) 
			A[k][i] = 0;

	return A_bidiag;

}

int main()
{
    
    int rows_A = 5;
    int cols_A = 5;
    
    double** A = matrix_read("csr_mat.txt", &rows_A, &cols_A);
    
    printf("rows and cols %d %d\n", rows_A, cols_A);
    printf("printing matrix\n");
    print_matrix(5, 5, A);

    double** gkb = bidiag(5, 5, A);
    
    printf("printing bidiagonal matrix\n");
    print_matrix(2, 5, gkb);

    printf("reflectors:\n");
    print_matrix(5, 5, A);
    
    return 0;
    
}