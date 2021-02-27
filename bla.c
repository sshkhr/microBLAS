#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bla.h"
#include "matvec_read.c"

// Question 6(a)
double inner(int length, double* vector_1, double* vector_2)
{
	int counter;
	double sum_of_products = 0;
	double vector_1_element, vector_2_element;

	for (counter = 0; counter < length; counter++)
	{
		vector_1_element = *(vector_1 + counter);
		vector_2_element = *(vector_2 + counter);
		sum_of_products += (vector_1_element * vector_2_element);
	}

	return sum_of_products;
}

// Question 6(b)
double norm(int length, double* vector)
{
	return p_norm(length, vector, 2);
}

// Question 6(c)
double normalize(int length, double* vector)
{
	double norm_value = norm(length, vector);

	divide_vector_by_number(length, vector, norm_value);

	return norm_value;
}

// Question 6(d)
double * normalize2(int length, double* vector)
{
	double norm_value = norm(length, vector);

	double * new_vector = (double *) malloc( length * sizeof(double) );

	for(int counter = 0; counter < length; counter++)
		*(new_vector + counter) = *(vector + counter);

	divide_vector_by_number(length, new_vector, norm_value);

	return new_vector;
}

// Question 6(e)
double * project(int length, double* vector_1, double* vector_2)
{

	double * projection_vector = (double *) malloc( length * sizeof(double) );

	projection_vector = normalize2(length, vector_2);

	double coefficient = inner(length, vector_1, vector_2)/norm(length, vector_2);

	multiply_vector_by_number(length, projection_vector, coefficient);

	return projection_vector;
}


// Question 6(g)
double * mat_vec_mult(int m, int n, double** matrix, double* vector)
{

	double* product_vector = (double*) malloc(m * sizeof(double) );

	int counter;
	
	for (counter = 0; counter < m; counter++)
	{
		double* matrix_row = *(matrix + counter);

		//printf("matrix row\n");
		//print_vector(4,matrix_row);

		*(product_vector + counter) = inner(n, matrix_row, vector);
	}

	return product_vector;
}


// Question 6(h)
double inner_col(int m, int i, double** matrix, double* vector)
{

	double* matrix_col = (double*) malloc( m * sizeof(double) );
	matrix_col = get_matrix_column(m, i, matrix);
	
	//printf("matrix col\n");
	//print_vector(3,matrix_col);

	return inner(m, matrix_col,vector);
}

// Question 6(i)
double inner_col2(int m, int col1, double** A, int col2, double** B)
{

	double* A_column = (double*) malloc( m * sizeof(double) );
	A_column = get_matrix_column(m, col1, A);

	//printf("matrix col\n");
	//print_vector(3,A_column);

	double* B_column = (double*) malloc( m * sizeof(double) );
	B_column = get_matrix_column(m, col2, B);

	//printf("matrix col\n");
	//print_vector(3,B_column);

	return inner(m, A_column, B_column);
}


int main()
{
	/*
	double arr[3] = {1,0,2};
	printf("norm %f\n",norm(3,arr));

	double new[3] = {2,0,3};
	printf("inner product %f\n",inner(3,arr,new));

	print_vector(3,new);
	printf("norm %f\n",normalize(3,new));
	print_vector(3,new);
	printf("norm %f\n",norm(3,new));

	double new_vec[3] = {2,0,3};
	
	print_vector(3,new_vec);
	double * new_vector = normalize2(3,new_vec);
	print_vector(3,new_vector);
	print_vector(3,new_vec);
	printf("norm %f\n",norm(3,new_vector));

	double vec_test[3] = {1,2,3};
	double vec_test_2[3] = {1,2,1};
	double * proj_vec = project(3,vec_test,vec_test_2);
	print_vector(3,proj_vec);

	int rows = 3;
	int cols = 4;
	double** mat = matrix_read("matrix1.txt", &rows, &cols);
	printf("printing matrix\n");
	print_matrix(3,4,mat);
	printf("inner col %f\n",inner_col(3,2,mat,vec_test_2));


	double** mat2 = matrix_read("matrix2.txt", &rows, &cols);
	printf("printing matrix\n");
	print_matrix(3,4,mat);
	printf("printing matrix\n");
	print_matrix(3,3,mat2);
	printf("inner col 2 %f\n",inner_col2(3,1,mat,2,mat2));
	
	double trial[4] = {1,2,3,4};
	double * mult_vector = mat_vec_mult(3, 4, mat, trial);
	print_vector(3,mult_vector);
	*/

	return 0;
}