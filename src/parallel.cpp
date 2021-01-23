#include <iostream>

#include "parallel.h"
#include "serial.h"
#include "util.h"
#include "omp.h"


/**
* @brief The paralle implementation of the **graph** algorithm.
* @param L: The sparse triangular matrix stored in the CSC format.
* @param b: The sparse vector stored in the CSC format.
* @param buffer: The buffer space in which the DFS operates.
* @param X: (Output parameter) The computed results.
* @return N/A.
*/
void lowTriSolveGraphPar(matrix_t<double> &L, matrix_t<double> &b, 
                      int *buffer, double *X) {

    int top_ptr = reachPar(L, b, buffer);

    // #pragma omp parallel for schedule(static)
    for (int i = 0; i < b.col_begins[b.N]; i++) {
        X[b.row_indices[i]] = b.values[i];
    }

    for (int i = top_ptr; i < L.N; i++) {
        int j = buffer[i];

        if (j < 0)
            continue;
        X[j] /= L.values[L.col_begins[j]];
        
        #pragma omp parallel for schedule(static)
        for (int val_ind = L.col_begins[j] + 1;
                 val_ind < L.col_begins[j + 1]; 
                 val_ind++) {
            X[L.row_indices[val_ind]] -= L.values[val_ind] * X[j];
        }
    }
}


/**
* @brief Perform depth-first search (DFS) on the dependency tree.
* @param start_node: The starting node of the search.
* @param start_ptr: The buffer pointer locates that locates its
*                   beginning position.
* @param buffer_upper: The upper portion of the buffer space 
*                      in which the recursion happens.
* @param buffer_lower: The lower portion of the buffer space 
*                      to which the nonzero-pattern goes.
* @param L: The sparse triangular matrix stored in the CSC format.
* @return The beginning of the buffer space.
*/
int dfsPar(int start, int top_ptr, int *buffer_upper, 
        int *buffer_lower, matrix_t<double> &L) {

    int buf_ptr = 0;
    bool is_finished = false;
    buffer_upper[buf_ptr] = start;
    
    while (buf_ptr >= 0) {
        is_finished = true;
        start = buffer_upper[buf_ptr]; 

        if (!IS_FLAGED(L.col_begins, start)) {
            FLAGE(L.col_begins, start);
            buffer_lower[buf_ptr] = (start < 0)? \
                    0 : RESTORE_SIGN(L.col_begins[start]);
        }

        int shift = (start < 0)? 0 : RESTORE_SIGN(L.col_begins[start+1]);
        for (int i = buffer_lower[buf_ptr]; i < shift; i++) {
            int val = L.row_indices[i];
            
            if (IS_FLAGED(L.col_begins, val))
                continue;
            
            is_finished = false;
            buffer_lower[buf_ptr] = i;
            buffer_upper[++buf_ptr] = val;
            break;
        }

        if (is_finished) {
            buf_ptr--;
            buffer_upper[--top_ptr] = start;
        }
    }
    return top_ptr;
}

/**
* @brief Compute the nonzero pattern of the matrix.
* @param L: The sparse triangular matrix stored in the CSC format.
* @param b: The sparse vector stored in the CSC format.
* @param buffer: The buffer space in which the DFS operates.
* @return The beginning of the buffer space.
*/
int reachPar(matrix_t<double> &L, matrix_t<double> &b, int *buffer) {
    int *buffer_upper = buffer + L.N;
    int top_ptr = L.N;
    // int t;

#pragma omp parallel
// #pragma omp single
{
//  lastprivate(top_ptr)
    #pragma omp for schedule(static)// lastprivate(t)
    for (int i = 0; i < b.col_begins[b.N]; i++) {
        #pragma omp critical
// #pragma omp task
        if (!IS_FLAGED(L.col_begins, b.row_indices[i])) {
            top_ptr = dfsPar(b.row_indices[i], top_ptr, buffer, buffer_upper, L);
            // top_ptr = t;
        }
    }
    // top_ptr = t;

} // END parallel
// #pragma omp barrier

    // #pragma omp parallel for schedule(static)
    for (int i = top_ptr; i < L.N; i++) {
        FLAGE(L.col_begins, buffer[i]);
    }

    return top_ptr;
}


/**
* @brief The parallel implementation of the **guarded** algorithm.
* @param dim: The matrix dimension.
* @param col_begins: The column pointer of L.
* @param row_indices: The row index of L.
* @param values: The values of L.
* @param bx: (I/O parameter) The right hand-side b at start_node 
*           and the solution x at the end.
* @return A EMPTY_INPUT(-1) if any of the inputs is empty 
*         or a EXIT_SUCCESS(1) when finish successfully.
*/
int lowTriSolveGuardedPar(int dim, int* col_begins,
                     int* row_indices, double* values, double* bx) {
    if (!col_begins || !row_indices || !bx)
        return EMPTY_INPUT;

    int j, p;
    #pragma omp parallel for schedule(static) private(p)
    for (j = 0; j < dim; j++) {
        if (bx[j] != 0) {
            bx[j] /= values[col_begins[j]];
            for (p = col_begins[j] + 1; p < col_begins[j + 1]; p++) {
                bx[row_indices[p]] -= values[p] * bx[j];
            }
        }
    }
    return EXIT_SUCCESS;
}


/**
* @brief The parallel implementation of the **naive** algorithm.
* @param dim: The matrix dimension.
* @param col_begins: The column pointer of L.
* @param row_indices: The row index of L.
* @param values: The values of L.
* @param bx: (I/O parameter) The right hand-side b at start_node 
*           and the solution x at the end.
* @return A EMPTY_INPUT(-1) if any of the inputs is empty 
*         or a EXIT_SUCCESS(1) when finish successfully.
*/
int lowTriSolveNaivePar(int dim, int* col_begins,
                     int* row_indices, double* values, double* bx) {
    if (!col_begins || !row_indices || !bx)
        return EMPTY_INPUT;

    for (int j = 0; j < dim; j++) {
        bx[j] /= values[col_begins[j]];

        #pragma omp parallel for schedule(static)
        for (int p = col_begins[j] + 1; p < col_begins[j + 1]; p++) {
            bx[row_indices[p]] -= values[p] * bx[j];
        }
    }
    return EXIT_SUCCESS;
}
