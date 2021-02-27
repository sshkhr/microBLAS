#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef BLA_H
#include "bla.c"
#endif

#ifndef READ
#include "matvec_read.c"
#endif

#include "qr_factorise.h"

double* scalar_product(int len, double s, double *v) {
	/* Compute the product of scalar s with vector v
	 * INPUT: length, scalar, vector
	 * RETURN: scaled vector */

	double *b;

	if ((b = (double *) calloc(len,sizeof(double))) == NULL) 
	{
		fprintf(stderr,"calloc failed\n");
		exit(1);
	}	
	int i;
	for (i=0; i<len; i++) b[i] = v[i]*s;
	return b;
}

double* vector_add(int len, double* v1, double* v2) {
	// Vector subtraction

	double *b;

	if ((b = (double *) calloc(len,sizeof(double))) == NULL) 
	{
		fprintf(stderr,"calloc failed\n");
		exit(1);
	}

	int i;
	for (i=0; i<len; i++) b[i] = v1[i] + v2[i];
	return b;
}

double* vector_sub(int len, double* v1, double* v2) {
	// Vector subtraction

	double *b;

	if ((b = (double *) calloc(len,sizeof(double))) == NULL) 
	{
		fprintf(stderr,"calloc failed\n");
		exit(1);
	}

	int i;
	for (i=0; i<len; i++) b[i] = v1[i] - v2[i];
	return b;
}


double** transpose(int row, int col, double **R) 
{
    int i, j;
    double temp;

    double** RT = (double **)malloc(col*sizeof(double*));
	for(int i = 0; i < col; i++)
	   RT[i] = (double *)malloc(row*sizeof(double));

    for (i = 0; i < col; i++) 
        for (j = 0; j < row; j++) 
        	RT[i][j] = R[j][i];

    return RT;
}

void qrhh1(int m, int n,  double** A)
{

	int i,j,p,q;

	for(i=0;i<n;i++)
	{
		// Allocate memory for A[i:m,i] as vector x and assign
		double* x = (double *)calloc(m-i,sizeof(double*));
		for(j=i;j<m;j++)
			x[j-i] = A[j][i];

		// Allocate memory for e1 and initialise
		double* e_one = (double *)calloc(m-i,sizeof(double*));
		e_one[0] = 1;

		//Get v_k = sign(x_1)norm(x)e_1 + x; v_k = v_k / ||v_k||;
		double* v_k = scalar_product(m-i, norm(m-i,x), e_one);
		v_k = vector_sub(m-i,x,v_k);
		//print_vector(m-i,v_k);		
		normalize(m,v_k);
		print_vector(m-i,v_k);		

		// Get reflector
		reflector = (double **)calloc(reflector,(m-i)*sizeof(double*));
		for(j = 0; j < (m-i); j++)
		   reflector[j] = (double *)calloc(reflector[j],(m-i)*sizeof(double));
		
		for(p=0;p<(m-i);p++)
			for(q=0;q<(m-i);q++)
				reflector[p][q] = 0;
		for(j = 0; j < (m-i); j++)
		   reflector[j][j] = 1;

		//printf("Allocated reflector\n");


		for(p=0;p<(m-i);p++)
			for(q=0;q<(m-i);q++)
				reflector[p][q] -= 2*v_k[p]*v_k[q];

		//printf("reflector\n");
		//print_matrix(m-i,m-i,reflector);

		// Get sub-matrix A[k:m][k:n]
		sub_matrix = (double **)calloc(sub_matrix,(m-i)*sizeof(double*));
		for(j = 0; j < (m-i); j++)
		   sub_matrix[j] = (double *)calloc(sub_matrix[j],(n-i)*sizeof(double));

		//printf("Allocated sub_matrix\n");

		for(p=0;p<(m-i);p++)
			for(q=0;q<(n-i);q++)
				sub_matrix[p][q] = A[p+i][q+i];
        
        //printf("sub matrix\n");
		//print_matrix(m-i,n-i,sub_matrix);

		//A[k:m][k:n] = A[k:m][k:n] - 2*v_k*(v_k*A[k:m][k:n])
		for(p=i;p<m;p++)
			for(q=i;q<n;q++)
				A[p][q] = inner_col(m-i, q, sub_matrix, reflector[p]);

		//print_matrix(m,n,A);

        // Free memory allocated to variables inside loop  
        free(x);
		free(e_one);
		free(v_k);   
		for(j = 0; j < (m-i); j++)
			free(reflector[j]);
		free(reflector);
		for(j = 0; j < (n-i); j++)
		   free(sub_matrix[j]);      
		free(sub_matrix);
		
	}

}

void qrhh2(int m, int n,  double** A, double* b)
{

	int i,j,p,q;

	double* Q_b = (double *)calloc(m,sizeof(double*));
	for(j=0;j<m;j++)
		Q_b[j] = b[j];		

	// take transpose s.t. column operations -> row operations
	//double** AT = transpose(m,n,A);

	for(i=0;i<n;i++)
	{
		// Allocate memory for A[i:m,i] as vector x and assign
		double* x = (double *)calloc(m-i,sizeof(double*));
		for(j=i;j<m;j++)
			x[j-i] = A[j][i];

		// Allocate memory for e1 and initialise
		double* e_one = (double *)calloc(m-i,sizeof(double*));
		e_one[0] = 1;

		//Get v_k = sign(x_1)norm(x)e_1 + x; v_k = v_k / ||v_k||;

		double* v_k = scalar_product(m-i, norm(m-i,x), e_one);
		v_k = vector_sub(m-i,x,v_k);
		//print_vector(m-i,v_k);		
		normalize(m,v_k);
		//print_vector(m-i,v_k);		

		// Get reflector
		double** reflector = (double **)calloc((m-i),sizeof(double*));
		for(j = 0; j < (m-i); j++)
		   reflector[j] = (double *)calloc((m-i),sizeof(double));
		for(j = 0; j < (m-i); j++)
		   reflector[j][j] = 1;

		for(p=0;p<(m-i);p++)
			for(q=0;q<(m-i);q++)
				reflector[p][q] -= 2*v_k[p]*v_k[q];

		//printf("reflector\n");
		//print_matrix(m-i,m-i,reflector);

		// Get sub-matrix A[k:m][k:n]
		double** sub_matrix = (double **)malloc((m-i)*sizeof(double*));
		for(j = 0; j < (m-i); j++)
		   sub_matrix[j] = (double *)malloc((n-i)*sizeof(double));

		for(p=0;p<(m-i);p++)
			for(q=0;q<(n-i);q++)
				sub_matrix[p][q] = A[p+i][q+i];
        
        //printf("sub matrix\n");
		//print_matrix(m-i,n-i,sub_matrix);

		//A[k:m][k:n] = A[k:m][k:n] - 2*v_k*(v_k*A[k:m][k:n])
		for(p=i;p<m;p++)
			for(q=i;q<n;q++)
				A[p][q] = inner_col(m-i, q, sub_matrix, reflector[p]);

		// Get Q*b
		Q_b = mat_vec_mult(m-i, n-i, reflector, Q_b);

		//print_matrix(m,n,A);

        // Free memory allocated to variables inside loop  
        free(x);
		free(e_one);
		free(v_k);   
		for(j = 0; j < (m-i); j++)
			free(reflector[j]);
		free(reflector);
		for(j = 0; j < (n-i); j++)
		   free(sub_matrix[j]);      
		free(sub_matrix);
		
	}

	for(j=0;j<m;j++)
		b[j] = Q_b[j];		

}
double** mgs1(int m, int n,  double** A)
{
	double** R;
	int i, j;
	double factor;
    
    // Allocate memory for R matrix
	R = (double **)calloc(n,sizeof(double*));
	for(int i = 0; i < n; i++) 
		R[i] = (double *)calloc(n,sizeof(double));
    
    // take transpose s.t. column operations -> row operations
	double** AT = transpose(m,n,A);

	for(i = 0; i<n; i++)
	{
		R[i][i] = norm(m, AT[i]);
		normalize(m, AT[i]);

		//printf("Row %d\n",i);
		//print_matrix(3,4,A);

		for(j = i+1; j<n; j++)
		{
			R[i][j] = inner(m,AT[i],AT[j]);
			double* sub = scalar_product(m,R[i][j],AT[i]);
			AT[j] = vector_sub(m,AT[j],sub);
		}

	}

	double** Q = transpose(n,m,AT);

	// Fill original matrix with Q matrix
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			A[i][j] = Q[i][j];

	//print_matrix(4,3,A);

	return R;
}

double** mgs2(int m, int n,  double** A, double* b)
{
	double** R;
	int i, j;
	double factor;
    
    // Allocate memory for R matrix
	R = (double **)calloc(n,sizeof(double*));
	for(i = 0; i < n; i++) 
		R[i] = (double *)calloc(n,sizeof(double));
    
    // Allocate memory for changed vector b_new
    double* b_new = (double *)calloc(m,sizeof(double));
    for(i = 0; i < m; i++) 
		b_new[i] = b[i];

    // take transpose s.t. column operations -> row operations
	double** AT = transpose(m,n,A);

	for(i = 0; i<n; i++)
	{
		R[i][i] = norm(m, AT[i]);
		normalize(m, AT[i]);

		//printf("Row %d\n",i);
		//print_matrix(3,4,A);

		for(j = i+1; j<n; j++)
		{
			R[i][j] = inner(m,AT[i],AT[j]);
			double* sub = scalar_product(m,R[i][j],AT[i]);
			AT[j] = vector_sub(m,AT[j],sub);
			b_new = vector_sub(m,b_new,sub);
		}

	}

	normalize(m, b);

	double** Q = transpose(n,m,AT);

	// Fill original matrix with Q matrix and vector b with b_new
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
			A[i][j] = Q[i][j];
		b[i] = b_new[i];
	}

	//print_matrix(4,3,A);
	//print_vector(4,b);
	return R;

}


/*
int main()
{
	
	int rows = 4;
	int cols = 3;
	double** mat = matrix_read("qr_test.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(4,3,mat);

	double** R = mgs1(4,3,mat);
	
    printf("printing Q matrix\n");
	print_matrix(4,3,mat);
	printf("printing R matrix\n");
	print_matrix(3,3,R);
	return 0;
	

	
	int rows = 3;
	int cols = 3;
	double** mat = matrix_read("qr_test_2.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(3,3,mat);

	double** R = mgs1(3,3,mat);
	
    printf("printing Q matrix\n");
	print_matrix(3,3,mat);
	printf("printing R matrix\n");
	print_matrix(3,3,R);
	return 0;
	

	
	int rows = 4;
	int cols = 3;
	double** mat = matrix_read("qr_test.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(4,3,mat);
	double vector_b[4] = {1,2,3,4};

	double** R = mgs2(4,3,mat,vector_b);
	
    printf("printing Q matrix\n");
	print_matrix(4,3,mat);
	printf("printing R matrix\n");
	print_matrix(3,3,R);
	printf("printing vector b\n");
	print_vector(4,vector_b);
	return 0;
	

	
	int rows = 4;
	int cols = 3;
	double** mat = matrix_read("qr_test.txt", &rows, &cols);
    printf("printing original matrix\n");
	print_matrix(4,3,mat);
	
	qrhh1(4,3,mat);
	
    printf("printing Q matrix\n");
	print_matrix(4,3,mat);
	return 0;
	
}
*/