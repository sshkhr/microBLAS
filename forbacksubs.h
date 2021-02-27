#ifndef FORBACKSUBS_H
#define FORBACKSUBS_H

double *backsubs(int n, double **R, double *b);
double *forsubs(int n, double **L, double *b);
double *backsubsT(int n, double **RT, double *b);
double *forsubsT(int n, double **LT, double *b);

#endif
