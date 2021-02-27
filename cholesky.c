#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "cholesky.h"

#ifndef BLA_H
#include "bla.c"
#endif

#ifndef READ
#include "matvec_read.c"
#endif

void cholesky(int m, double** A)
{
	double squared_sum = 0;
	double sum = 0;

	for(int i=0;i<m;i++) //iterate over columns
	{
		if(i) // from second column onwards
			squared_sum = inner_col2(i, i, A, i, A);

		A[i][i] = sqrt(A[i][i] - squared_sum);

		for(int j=i+1;j<m;j++)
		{
			if(i) // from second column onwards
				sum = inner_col2(i, j, A, i, A);

			// upper triangular term
			A[i][j] = (A[i][j] - sum)/A[i][i];

			// lower triangular term
			A[j][i] = 0;
		}
	}
}

/*
int main()
{
	
	int rows = 5;
	int cols = 5;
	double** mat = matrix_read("cholesky_test_new_5.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(5,5,mat);

	cholesky(5,mat);
	
    printf("printing R matrix\n");
	print_matrix(5,5,mat);
	return 0;
}
*/