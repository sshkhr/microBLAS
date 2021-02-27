#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matvec_read.h"

double **matrix_read(char *filename, int *rows, int *cols) {
	/* Read a text file that defines a matrix.
	 * Each row should have the same number of entries as the first row. 
	 * Additional entries are ignored; missing entries are treated as zero.
	 * Entries in a row are separated by whitespace or commas.
	 *
	 * Arguments: 
	 *		input file name,
	 *		pointer to int to store number of rows
	 *		pointer to int to store number of columns
	 *
	 * Returns:  
	 *		NULL pointer on failure
	 *		Pointer to pointer to double otherwise.  (Array of pointers to rows.)
	 */

	int m,n,rowlen;
	double **matrix;
	char b;
	char *line,*tok;
	char delimiters[5] = ", \t\n";
	FILE *infile;

	if ((infile = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Failed to open file %s.\n",filename);
		return NULL;
	}
	/* First read characters until the first newline character is encountered. */
	rowlen = 0;
	b = fgetc(infile);
	while (b != EOF && b != '\n') {
		rowlen++;
		b = fgetc(infile);
	}
	rowlen += 2;
	rewind(infile);
	/* Allocate space for, and read the first line into a string. */
	line = (char *) malloc(rowlen*sizeof(char));
	fgets(line,rowlen,infile);
	/* Count the number of entries in this line. */
	tok = strtok(line,delimiters);
	*cols = 0;
	while (tok != NULL) {
		(*cols)++;
		tok = strtok(NULL,delimiters);
	}
	free(line);
	/* Count the remaining lines, keeping track of the longest line. */
	*rows = 1;
	n = 0;
	b = fgetc(infile);
	while (b != EOF) {
		if (b == '\n') {
			(*rows)++;
			if (n > rowlen) rowlen = n;
		}
		else
			n++;
		b = fgetc(infile);
	}
	/* Allocate space for longest line, and for row pointers. */
	if (((line = (char *) malloc((rowlen+2)*sizeof(char))) == NULL) || 
	    ((matrix = (double **) malloc(*rows*sizeof(double *))) == NULL)) {
		fprintf(stderr,"malloc failed\n");
		fclose(infile);
		return NULL;
	}
	/* Rewind file and read lines. */
	rewind(infile);
	/* Allocate space for all rows. */
	for (m=0; m<*rows; m++) {
		if ((matrix[m] = (double *) calloc(*cols,sizeof(double))) == NULL) {
			fprintf(stderr,"calloc failed\n");
			fclose(infile);
			/* error getting row m; free previous rows */
			for (n=m-1; n>=0; n--) free(matrix[n]);
			free(matrix);
			free(line);
			return NULL;
		}
	}
	/* Read each line of file and record values. */
	for (m=0; m<*rows; m++) {
		/* Get the line from the file. */
		if (fgets(line,rowlen,infile) == NULL) {
			fprintf(stderr,"Error reading input file line %d.\n",m+1);
			m = *rows*10;  /* to indicate error condition */
		}
		else {
			/* Break the line by white space. */
			tok = strtok(line,delimiters);
			for (n=0; n<*cols; n++) {
				if (tok == NULL) break;  /* missing columns treated as zeros */
				if (sscanf(tok,"%lf",matrix[m]+n) != 1) {
					fprintf(stderr,"Error on line %d column %d of input file %s:%s%s.\n",
							m+1,n+1,filename," Trying to read: ",tok);
					m = *rows*10;  /* to indicate error condition */
				}
				else tok = strtok(NULL,delimiters);
			}
		}
	}
	if (m == *rows*10) {   /* error condition occurred */
		/* free space */
		for (m=0; m<*rows; m++) free(matrix[m]);
		free(matrix);
		matrix = NULL;
	}
	free(line);
	fclose(infile);
	return matrix;
}

double *vector_read(char *filename, int *entries) {
	/* Read a text file that defines a vector.
	 * The file is read row-by-row, taking as many entries as are available on 
	 * each line and putting them all together in one vector, white space and commas
	 * are ignored.
	 *
	 * Arguments: 
	 *		input file name,
	 *		pointer to int to store number of entries
	 *
	 * Returns:  
	 *		NULL pointer on failure
	 *		Pointer to double otherwise.
	 */

	int n;
	double *vector,*newvector;
	FILE *infile;

	if ((infile = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Failed to open file %s.\n",filename);
		return NULL;
	}
	/* Allocate space starting with 50 units, then either doubling or adding 1000 each time more is
	 * needed. */
	*entries = 50;
	if ((vector = (double *) malloc(*entries*sizeof(double))) == NULL) {
		fprintf(stderr,"Error allocating initial space for vector.\n");
		return NULL;
	}
	n = 0;
	while (fscanf(infile,"%lf%*[, \t\r]",vector+n) == 1) {
		n++;
		if (n == *entries) {
			if (*entries < 1000)
				*entries *= 2;
			else
				*entries += 1000;
			if ((newvector = (double *) realloc(vector, *entries*sizeof(double))) == NULL) {
				fprintf(stderr,"Error increasing space to %d for vector.\n",*entries);
				free(vector);
				return NULL;
			}
			vector = newvector;
		}
	}
	if (*entries != n) {
		if ((newvector = (double *) realloc(vector, n*sizeof(double))) == NULL) 
			fprintf(stderr,"Warning: failed to reduce vector to size %d.\n",n);
		else vector = newvector;
		*entries = n;
	}
	fclose(infile);
	return vector;
}

