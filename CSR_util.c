#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CSR_util.h"
#include "matvec_read.c"

CSR* full2CSR(int rows, int cols, double** matrix)
{
	CSR* csr_matrix = malloc(sizeof(CSR));

	csr_matrix->m = (int *) malloc (sizeof *csr_matrix->m);
	csr_matrix->n = (int *) malloc (sizeof *csr_matrix->n);
	
	*csr_matrix->m = rows;
	*csr_matrix->n = cols;

	int row_counter;
	int col_counter;
	double value;
	double* row_address;

	// count number of non-zeros
	int non_zero = 0;

	for (row_counter = 0; row_counter < *(csr_matrix->m); row_counter++)
	{
		row_address = *(matrix + row_counter);

		for (col_counter = 0; col_counter < *(csr_matrix->n); col_counter++)
		{
			value = *(row_address + col_counter);

			if(value != 0)
				non_zero++;
		}
	}


	// Memory allocation for V, J, I
	csr_matrix->V = (double *) malloc (non_zero * sizeof (*csr_matrix->V));
	csr_matrix->J = (double *) malloc (non_zero * sizeof (*csr_matrix->J));
	csr_matrix->I = (double *) malloc ((*(csr_matrix->m)+1) * sizeof (*csr_matrix->I));

	double* original_v = (double *) malloc (sizeof(double));
	double* original_j = (double *) malloc (sizeof(double));
	double* original_i = (double *) malloc (sizeof(double));

	// Save the starting pointers for V, J, I
	original_v = csr_matrix->V;
	original_j = csr_matrix->J;
	original_i = csr_matrix->I;

	non_zero = 0;

	for (row_counter = 0; row_counter < *(csr_matrix->m); row_counter++)
	{
		// Store number of non-zeros encountered till this row
		*(csr_matrix->I) = non_zero;

		row_address = *(matrix + row_counter);

		for (col_counter = 0; col_counter < *(csr_matrix->n); col_counter++)
		{
			
			value = *(row_address + col_counter);

			if(value != 0)
			{
				non_zero++;
				
				// Save value and column in V,J vectors and increment pointers

				*csr_matrix->V = value;
				*csr_matrix->J = col_counter;

				csr_matrix->V += 1;
				csr_matrix->J += 1;
				
			}
		}
		
		// Increment pointer for i vector
		csr_matrix->I += 1;
	}

	// Store total number of non zeros as last element of I
	*(csr_matrix->I) = non_zero;

	// Restore the starting pointers for V, J, I
	csr_matrix->V = original_v;
	csr_matrix->J = original_j;
	csr_matrix->I = original_i;

	return csr_matrix;
}

double* CSRmult(CSR* csr_matrix, double vector[])
{
	int m = *(csr_matrix->m);
	int n = *(csr_matrix->n);
	
	double* product_vector = (double*) malloc(m * sizeof(double) );

	int row_counter, col_counter, column_index;
	double matrix_element, vector_element;

	// Initialise the product vector values as 0
	for (row_counter = 0; row_counter < m; row_counter++)
		product_vector[row_counter] = 0;
	
	// For matrix vector multiplication as row of matrix * vector
	// Row can be kept track of by using the I vector in CSR

	for (row_counter = 0; row_counter < m; row_counter++)
	{
		// This loop tracks the number of non-zero elements in this row
		for(col_counter = (int)*(csr_matrix->I + row_counter); col_counter < (int)*(csr_matrix->I + row_counter + 1); col_counter++)
		{
			// The multiplication has following steps:

			// *(csr_matrix->V + col_counter) gets value of non-zero element A[row_counter][column_counter]
			matrix_element = *(csr_matrix->V + col_counter);
			
			// Get the column_index of non-zero element and corresponding vector element
			column_index = (int) *(csr_matrix->J + col_counter);
			vector_element = vector[column_index];

			//printf("matrix element %f * vector element %f \n",matrix_element, vector_element);

			// The prduct of these two numbers is added to product vector
			// This process gets looped over all non-zero elements in row to get row*vector value

			product_vector[row_counter] += (matrix_element * vector_element);
		}
	}

	return product_vector;
}

int main()
{
	return 0;
    
    /*
	int rows = 5;
	int cols = 5;
	double** mat = matrix_read("csr_mat.txt", &rows, &cols);
	printf("rows and cols %d %d\n", rows, cols);
	printf("printing matrix\n");
	print_matrix(5, 5, mat);

	CSR* out = full2CSR(5, 5, mat);
	printCSR(out);

	double vec_test[5] = {1,2,3,4,5};
	double * product = CSRmult(out, vec_test);
	print_vector(5,product);
	*/
}