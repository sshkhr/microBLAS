#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "eigenvalues.h"

#ifndef BLA_H
#include "bla.c"
#endif

#ifndef QRR_H
#include "qrhh3.c"
#endif

#ifndef READ
#include "matvec_read.c"
#endif

void upphess(int m, double **A)
{
	int i,j,k;
	double v[m-1],normv,sx1,vw;

	for(k=0;k<m-1;k++) //iterate over columns
	{

		// Vector that x is projected on to
		for(i=0;i<m-k-1;i++)
			v[i] = A[i+k+1][i];

		// Squared norm of column i below sub-diagonal
		normv = inner(m-k-2,v+1,v+1);

		//if norm zero don't need to do anything for the column
		if(normv)
		{
			//Get the normalised vector v
			normv = sqrt(normv + v[0]*v[0]);
			sx1 = (v[0] >= 0.0) ? normv : -normv;
			v[0] += sx1;
			normv = normalize(m-k-1,v);
		
			A[k+1][k] = -sx1;

			// zero out current column below sub diagonal
			for(i=k+2; i<m; i++)
				A[i][k] = 0;

			// Apply reflector on left
			for(j=k+1;j<m;j++)//columns to the right
			{
				vw = inner_col(m-k-1,j,&(A[k+1]),v);
				for(i=0;i<m-k-1;i++) // rows below subdiagonal
					A[i+k+1][j] -= 2*v[i]*vw;
			}

			// Apply reflector on right
			for(i=0;i<m;i++) // all rows
			{
				vw = inner(m-k-1,&(A[i][k+1]),v);
				for(j=0;j<m-k-1;j++) // columns to the right
					A[i][j+k+1] -= 2*v[j]*vw;
			}			
		}

	}

}

int qralg(int n, double **A)
{
	double tolerance = 1.0E-14;
	double a,b,c,d,sign,delta,mu;

	int iteration = 0;
	int i;

	while(abs(A[n-1][n-2]) > tolerance)
	{

		// Get the shift term mu
		a = A[n-2][n-2];
		b = A[n-2][n-1];
		c = A[n-1][n-2]; //b and c are same
		d = A[n-1][n-1];

		delta = (a-d)/2;
		sign = (delta > 0)? 1 : -1;
		mu = d - sign*pow(b,2)/(abs(delta) + sqrt(pow(delta,2) + pow(b,2)));

		// Get A-mu.I
		for(i=0;i<n;i++)
			A[i][i] -= mu;

		// Get QR factorisation of shifted A
		qrhh3(n, A);

		// Update A
		for(i=0;i<n;i++)
			A[i][i] += mu;

		iteration++;	
	}

	printf("QR complete\n");
	print_matrix(n,n,A);

	return iteration;
	
}

int* eigval(int m, double **A)
{
	int* iterations;

	if ((iterations = (int *) calloc(m,sizeof(int))) == NULL) 
	{
		fprintf(stderr,"calloc failed\n");
		exit(1);
	}

	upphess(m,A);

	for(int i=0;i<m;i++)
		iterations[i] = qralg(m-i, A);

	return iterations;
}

/*
int main()
{
	
	int rows = 5;
	int cols = 5;
	double** mat = matrix_read("eig_test.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(5,5,mat);

	int* P = eigval(5,mat);
	
    printf("printing diagnal matrix\n");
	print_matrix(5,5,mat);
	printf("printing iterations vector\n");
	
	printf("vector p: \t");
		for(int q=0;q<5;q++)
			printf("%d\t",P[q]);
	printf("\n");

	return 0;

	int rows = 10;
	int cols = 10;
	double** mat = matrix_read("eigval_test_matrix.txt", &rows, &cols);
    
	int* P = eigval(10,mat);
	
    printf("printing diagnal matrix\n");
	print_matrix(10,10,mat);
	printf("printing iterations vector\n");
	
	printf("vector p: \t");
		for(int q=0;q<10;q++)
			printf("%d\t",P[q]);
	printf("\n");

	return 0;

}
*/