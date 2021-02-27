#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bla.h"
#include "CSR_util.h"
#include "arnoldi.h"
#include "matvec_read.h"

#include "helper.c"

#define MAX(a, b) ((a) >= (b) ? (a) : (b))

int arnoldi(struct CSR *A, int n, double ***QT, double ***HT, double *b)
{
	static int iterations = 0;
	static int rows_allocated = 0;
	int rows_to_allocate;

	int rows_A = A->m;
	int cols_A = A->n;

	int i,j,k;

	double *v;

	// use b to initialize QT in first iteration
	if(iterations == 0)
	{
		rows_allocated += MAX(3*n,10);

		*QT = calloc(rows_allocated,sizeof(double *));
		for (i=0; i<rows_allocated; i++)
		    *(*QT+i) = calloc(rows_A,sizeof(double));

		*HT = calloc(rows_allocated,sizeof(double *));
		for (i=0; i<rows_allocated; i++)
		    *(*HT+i) = calloc(i+2,sizeof(double));

		normalize(rows_A,b);

		for (i=0; i<rows_A; i++) 
			*(*(*QT)+i) = b[i];
	}

	// perform arnoldi n times
	for(i = 0; i < n; iterations++)
	{
		if(iterations == rows_A-1)
			return i;

		// Check if more rows need to be allocated
		if(iterations == rows_allocated)
		{
            rows_to_allocate = MAX(3*n,10);
            rows_allocated += rows_to_allocate;

            /* Reallocate QT */
			*QT = realloc(*QT,rows_allocated*sizeof(double *));
			for (j=rows_allocated-rows_to_allocate; j<rows_allocated; j++)
			    *(*QT+j) = calloc(rows_A,sizeof(double));

			/* Reallocate HT */
			*HT = realloc(*HT,rows_allocated*sizeof(double *));
			for (j=rows_allocated-rows_to_allocate; j<rows_allocated; j++)
			    *(*HT+j) = calloc(j+2,sizeof(double));
		}

		// v = A*Qn
		v = CSRmult(A, *(*QT+iterations));

		for(j=0;j<iterations;j++)
		{
			// H[j][i] = Q[j].v
			*(*(*HT+iterations)+j) = inner(rows_A,*(*QT+iterations),v);

			// v = v - H[j][i]*Q[j] 
			for(k=0;k<rows_A;k++)
				v[k] = v[k] - *(*(*HT+iterations)+k) * *(*(*QT+iterations)+k);
		}

		// H[i+1][i] = ||v||
		*(*(*HT+iterations)+iterations+1) = normalize(rows_A,v);

		// Arnoldi breaks down. Return #iterations run
		if(*(*(*HT+iterations)+iterations+1) < (rows_A*1E-15)) 
			return i;

		// Q[i+1] = v_normalized
		for(k=0;k<rows_A;k++)
			*(*(*QT+iterations+1)+k) = v[k];

		i++;
	}

	return n;
}

int main()
{
	int rows_A = 5;
    int cols_A = 5;
    
    double** A = matrix_read("csr_mat.txt", &rows_A, &cols_A);
    
    //printf("rows and cols %d %d\n", rows_A, cols_A);
    //printf("printing matrix\n");
    //print_matrix(5, 5, A);

    struct CSR* csr_A = full2CSR(5,5,A);

    //printCSR(csr_A);

	double **QT; 
	double **HT;

	double b[5] = {1,2,3,4,5}; 

	//print_vector(5,b);

	int first = arnoldi(csr_A, 1, &QT, &HT, b);
    
    printf("Run for %d\n", first);

    printf("QT\n");
    print_matrix(2, 5, QT);
    
	int second = arnoldi(csr_A, 2, &QT, &HT, NULL);
    
    printf("Run for %d\n", second);


    printf("QT\n");
    print_matrix(4, 5, QT);
    
    int third = arnoldi(csr_A, 4, &QT, &HT, NULL);
    
    printf("Run for %d\n", third);

    
    printf("QT\n");
    print_matrix(5, 5, QT);
    
    return 0;    
}