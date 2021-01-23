#ifndef _PARALLEL_H
#define _PARALLEL_H

#include <iostream>
#include "serial.h"

void lowTriSolveGraphPar(matrix_t<double> &L, matrix_t<double> &b, int *xi, double *x);
int  lowTriSolveNaivePar(int dim, int *col_begins, int *row_indices, double *values, double *bx);
int  lowTriSolveGuardedPar(int dim, int *col_begins, int *row_indices, double *values, double *bx);

int reachPar(matrix_t<double> &L, matrix_t<double> &b, int *buffer);
int dfsPar(int start, int top_ptr, int *buffer_upper, 
           int *buffer_lower, matrix_t<double> &L);

#endif // _PARALLEL_H