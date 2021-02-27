#include <stdlib.h>
#include <stdio.h>

#include "CSRmat_mult.h"
#include "CSR_util.h"
#include "matvec_read.h"

#include "helper.c"

struct CSR *CSRmat_mult(struct CSR *A, struct CSR *B)
{
    struct CSR *prod;

    /* Get space for the CSR structure. */
    if ((prod = (struct CSR *) malloc(sizeof(struct CSR))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    }

    int m = A->m;
    int n = B->n;

    prod->m = m;
    prod->n = n;

    /* Get space for I in the CSR structure. */
    if ((prod->I = (int *) calloc(m+1,sizeof(int))) == NULL) {
        fprintf(stderr,"calloc failed\n");
        exit(1);
    }

    /* Get space for V and J in the CSR structure. */
    if ((prod->V = (double *) malloc((prod->m)*sizeof(double))) == NULL ||
            (prod->J = (int *) malloc((prod->m)*sizeof(int))) == NULL) 
    {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    }

    int non_zero_A = A->I[m];
    int non_zero_B = B->I[m];
    double prod_i_j, a_element, b_element;

    int a_row_start, a_row_end, b_row_start, b_row_end, a_col_iterator, b_col_iterator;
    int a_row, b_row, a_col, b_col, prod_row, prod_col;

    int non_zero_row, index;

    for(prod_row = 0; prod_row < prod->m; prod_row++)
    {
        non_zero_row = 0; // flag to check for non-zero elements in row

        for(prod_col = 0; prod_col < prod->n; prod_col++)
        {
            // calculate element prod[prod_row][product_col]
            // prod_i_j = A[i][:]*B[:][j] 
            // prod_i_j += A[i][k]*B[k][j]

            prod_i_j = 0;

            a_row = prod_row;

            a_row_start = A->I[a_row];
            a_row_end = A->I[a_row+1];

            if(a_row_start == a_row_end) // row of zeros in A
                continue; // denotes prod_i_j = 0

            // iterate over all the columns in row i of A
            for(a_col_iterator=a_row_start;a_col_iterator<a_row_end;a_col_iterator++)
            {
                a_col = A->J[a_col_iterator]; // get column number of ith row element A[i][k]
                a_element = A->V[a_col_iterator]; // get A[i][k]

                // Iterate over all elements in B
                for(b_col_iterator=0;b_col_iterator<non_zero_B;b_col_iterator++)
                {
                    b_col = B->J[b_col_iterator];

                    // if column of B not same as product column continue
                    if(b_col != prod_col)
                    {
                        continue;
                    }
                    else
                    {
                        // Have A[i][k] and B[:][j]
                        // check if row of B is same as column of A (k)

                        b_row_start = B->I[a_col];
                        b_row_end = B->I[a_col+1];

                        // if the current element belongs to this row
                        if(b_col_iterator>=b_row_start && b_col_iterator<b_row_end)
                        {
                            b_element = B->V[b_col_iterator];

                            // Multiply A[i][k] and B[k][j]
                            prod_i_j += a_element*b_element;
                        }

                    }
                }

            }
            
            if(prod_i_j == 0)
                continue;
            else
            {
                non_zero_row += 1; // increment # non-zero elements in this row

                // get index in CSR format
                index = prod->I[prod_row] + non_zero_row - 1;

                /* Reallocate space for V and J in the CSR structure. */
                if (((prod->V = realloc(prod->V,(index+1)*sizeof(double))) == NULL) ||
                        ((prod->J = realloc(prod->J,(index+1)*sizeof(int))) == NULL)) 
                {
                    fprintf(stderr,"realloc failed\n");
                    exit(1);
                }

                prod->J[index] = prod_col;
                prod->V[index] = prod_i_j;
                
            }
        }

        // start point of next row is starting point of previous + # non-zero elements in this row
        prod->I[prod_row+1] = prod->I[prod_row] + non_zero_row;
    }

    return prod;

}


int main()
{
    
    int rows_A = 5;
    int cols_A = 5;
    double** A = matrix_read("csr_mat.txt", &rows_A, &cols_A);
    //printf("rows and cols %d %d\n", rows, cols);
    //printf("printing matrix\n");
    //print_matrix(5, 5, mat);

    struct CSR* out_A = full2CSR(5, 5, A);
    printCSR(out_A);

    int rows_B = 5;
    int cols_B = 2;
    double** B = matrix_read("csr_mat_2.txt", &rows_B, &cols_B);
    //printf("rows and cols %d %d\n", rows, cols);
    //printf("printing matrix\n");
    //print_matrix(5, 5, mat);

    struct CSR* out_B = full2CSR(5, 2, B);
    printCSR(out_B);

    struct CSR* mult = CSRmat_mult(out_A, out_B);
    printCSR(mult);

    return 0;
    
}