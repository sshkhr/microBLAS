#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double inner(int length, double* vector_1, double* vector_2);
double norm(int length, double* vector);
double normalize(int length, double* vector);
double * normalize2(int length, double* vector);
double * project(int length, double* vector_1, double* vector_2);
double * mat_vec_mult(int m, int n, double** matrix, double* vector);
double inner_col(int m, int i, double** matrix, double* vector);
double inner_col2(int m, int col1, double** A, int col2, double** B);


// This helper function prints out a vector for debugging
void print_vector(int length, double* vector)
{
	int counter;
	
	for (counter = 0; counter < length; counter++)
		printf("%f\t",*(vector + counter));

	printf("\n");
}

// This helper function prints out a matrix for debugging
void print_matrix(int rows, int cols, double** matrix)
{
	printf("%d, %d\n", rows, cols);
	int row_counter;
	int col_counter;
	
	for (row_counter = 0; row_counter < rows; row_counter++)
	{
		printf("row counter %d \n", row_counter);
		double* row_address = *(matrix+row_counter);
		printf("row_address %p \n",(double *)row_address);
		for (col_counter = 0; col_counter < cols; col_counter++)
		{

			printf("%f\t",*(row_address+col_counter));
		}
		printf("\n");
	}
	printf("\n");
}

// This helper function computes and returns the p-norm of a vector
double p_norm(int length, double* vector, double p)
{
	int counter;
	double sum_of_powers = 0;
	double vector_element;

	for (counter = 0; counter < length; counter++)
	{
		vector_element = *(vector + counter);
		sum_of_powers += pow(vector_element, p);
	}

	double norm = pow(sum_of_powers, (1/p));

	return norm;
}

// This helper function divides a vector of doubles of given length by a scalar in-place
void divide_vector_by_number(int length, double* vector, double constant)
{
	int counter;
	double vector_element;

	for(counter = 0; counter < length; counter++)
	{
		if(constant)
			*(vector + counter) = *(vector + counter)/constant;
		else
			*(vector + counter) = *(vector + counter);
	}
}

// This helper function multiplies a vector of doubles of given length by a scalar in-place
void multiply_vector_by_number(int length, double* vector, double constant)
{
	int counter;
	double vector_element;

	for(counter = 0; counter < length; counter++)
			*(vector + counter) = (*(vector + counter))*constant;
}


// This helper function extracts a particular column from matrix stored as pointer to pointer to oduble
double* get_matrix_column(int rows, int col, double** matrix)
{
	int counter;
	double* matrix_col = (double*) malloc( rows * sizeof(double) );

	for (counter = 0; counter < rows; counter++)
	{
		double* row_address = *(matrix+counter);
		*(matrix_col + counter) = *(row_address + col);
	}
		
	return matrix_col;
}