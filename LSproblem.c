#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifndef BLA_H
#include "bla.c"
#endif

#ifndef READ
#include "matvec_read.c"
#endif


#include "LSproblem.h"
#include "qr_factorise.c"
#include "forbacksubs.c"

#define pi 3.141592654

void AB_construct(int n, double*** A, double** b)
{
	int i;
	double x;

	*A = (double **)calloc(n,sizeof(double*));

	*b = (double *)calloc(n,sizeof(double)); 

	for(i = 0; i < n; i++)
		*(*A+i) = (double *)calloc(4,sizeof(double));
	
	for(i = 0; i < n; i++)
	{

		x = ((double)(i+1))/n;

		*(*(*A+i)+0) = 1;
		*(*(*A+i)+1) = pow(x,1);
		*(*(*A+i)+2) = pow(x,2);
		*(*(*A+i)+3) = pow(x,3);

		*(*b+i) = cos(3*pi*x);		
	}

	//print_matrix(5,4,*A);

}

double LS_solve(int n, double** c)
{
	double** A;
	double* b;

	AB_construct(n,&A,&b);

	//print_matrix(n,4,A);
	//print_vector(n,b);

	double** R = mgs2(n,4,A,b);

	//print_matrix(4,4,R);
	//print_vector(4,b);
	
	*c = (double *)calloc(4,sizeof(double));

	*c = backsubs(4, R, b);

	double x,residual;
	double E = 0;
	double x_row[4];

	for(int i = 0; i < n; i++)
	{
		x = ((double)(i+1))/n;

		x_row[0]=1;x_row[1]=pow(x,1);x_row[2]=pow(x,2);x_row[3]=pow(x,3);

		residual = cos(3*pi*x) - inner(4, *c, x_row);

		E += pow(residual,2);
	}

	E /= n;
	E = sqrt(E);
	return E;
}

/*
int main()
{
	
	double** A;
	double* b;

	A = (double **)calloc(5,sizeof(double*));
	for(int i = 0; i < 5; i++)
		A[i] = (double *)calloc(4,sizeof(double));

	b = (double *)calloc(5,sizeof(double));

	AB_construct(4,&A,&b);

	print_matrix(4,4,A);

	print_vector(4,b);
	
	double* c;
	double E;

	double E_prev = LS_solve(4,&c);
	E_prev = round( E_prev * 1000000.0 ) / 1000000.0;

	int i=5;

	while(i)
	{
		E =  LS_solve(i,&c);
		E = round( E * 1000000.0 ) / 1000000.0;

		if(E_prev == E)
		{
			printf("n = %d E = %f\n",i,E);
			printf("c =");
			print_vector(4,c);
			break;
		}
		else
			E_prev = E;

		i++;
	}

	return 0;
}
*/