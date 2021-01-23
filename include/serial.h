#ifndef _SOLVE_H
#define _SOLVE_H

#include <iostream>

#define NEGATE(i) (-(i)-0xffff) 
#define IS_FLAGED(arr,i) (arr[i] < 0) 
#define RESTORE_SIGN(i) (i < 0) ? NEGATE(i) : i
#define FLAGE(arr,i) {arr[i] = NEGATE(arr[i]);}

template<class T>  struct matrix_t {
    int N;      // #Columns
    int M;      // #Rows
    int nnz;    // #Nonezeros
    
    /* CSC format */
    T   *values;
    int *col_begins;
    int *row_indices;
};

typedef struct matrix_t<double> mtx_dbl;

void   lowTriSolveGraph(matrix_t<double> &L, matrix_t<double> &b, int *xi, double *x);
size_t lowTriSolveNaive(int dim, int *col_begins, int *row_indices, double *values, double *bx);
size_t lowTriSolveGuarded(int dim, int *col_begins, int *row_indices, double *values, double *bx);

int reach(matrix_t<double> &L, matrix_t<double> &b, int *buffer);
int dfs(int start, int top_ptr, int *buffer_upper, 
        int *buffer_lower, matrix_t<double> &L);

#endif // _SOLVE_H