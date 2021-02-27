#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "palu.h"

#ifndef BLA_H
#include "bla.c"
#endif

#ifndef READ
#include "matvec_read.c"
#endif

void swap_row(double** A, int i, int j, int m) 
{ 
    int k;
    double temp;
  
    for (k=0; k<m; k++) 
    { 
        temp = A[i][k]; 
        A[i][k] = A[j][k]; 
        A[j][k] = temp; 
    }
} 

int* palu(int m, double** A)
{
	int i,j,k, max_row, *p, temp;
	double max_element, ratio;

	if ((p = (int *) calloc(m,sizeof(int))) == NULL) 
	{
		fprintf(stderr,"calloc failed\n");
		exit(1);
	}

	for(int q=0;q<m;q++)
		p[q] = q;
		

	
	for(j=0;j<m;j++) //iterate over columns
	{
		max_row = j;
		max_element = A[max_row][j];

		for(i=j+1;i<m;i++) //iterate over all rows below column number
		{
			if(abs(A[i][j]) > max_element)
				{
					max_row = i;
					max_element = abs(A[i][j]);
				}
		}

		if(!A[max_row][j])
		{
			// Singular
		}

		if(max_row != j)
		{
			printf("Swapping rows %d and %d\n",max_row,j);
			swap_row(A, max_row, j, m);


			// Change the permutation sequence in p
			temp = p[j];
			p[j] = p[max_row];
			p[max_row] = temp;
			
		}

		print_matrix(m,m,A);
		
		printf("vector p: \t");
		for(int q=0;q<m;q++)
			printf("%d\t",p[q]);
		printf("\n");

		// iterate over all rows 
		// divide by ratio, subtract with pivot row and get zeros
		for(i=j+1;i<m;i++) 
		{
			ratio = A[i][j]/A[j][j];

			for(k=j+1;k<m;k++)
				A[i][k] -= ratio*A[j][k];

			A[i][j] = ratio; // get element of lower traingular L 
		}
	}

	
	return p;
}

/*
int main()
{
	
	int rows = 4;
	int cols = 4;
	double** mat = matrix_read("palu_test.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(4,4,mat);

	int* P = palu(4,mat);
	
    printf("printing L-U matrix\n");
	print_matrix(4,4,mat);
	printf("printing P vector\n");
	
	printf("vector p: \t");
		for(int q=0;q<4;q++)
			printf("%d\t",P[q]);
	printf("\n");

	return 0;
}
*/