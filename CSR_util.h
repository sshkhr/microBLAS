#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bla.h"

typedef struct 
{
	int* m;
	int* n;

	double* V;
	double* J;
	double* I;
} CSR;

CSR* full2CSR(int m, int n, double** matrix);
double* CSRmult(CSR* csr_matrix, double vector[]);


// this helper function prints out a matrix in CSR format for debugging
void printCSR(CSR* csr_matrix)
{
	printf("rows %d\n",*(csr_matrix->m));
	printf("cols %d\n",*(csr_matrix->n));

	int non_zeros = *(csr_matrix->I + *(csr_matrix->m));
	printf("Non zeros %d\n",non_zeros);

	printf("V vector\n");
	print_vector(non_zeros,csr_matrix->V);

	printf("J vector\n");
	print_vector(non_zeros,csr_matrix->J);

	printf("I vector\n");
	print_vector(*(csr_matrix->m)+1,csr_matrix->I);
}