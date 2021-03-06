#include <iostream>

#include "serial.h"
#include "util.h"
#include "omp.h"


/**
* @brief The **graph** algorithm for solving the triangular system Lx=b.
* @param L: The sparse triangular matrix stored in the CSC format.
* @param b: The sparse vector stored in the CSC format.
* @param buffer: The buffer space in which the DFS operates.
* @param X: (Output parameter) The computed results.
* @return N/A.
*/
void lowTriSolveGraph(matrix_t<double> &L, matrix_t<double> &b, 
                      int *buffer, double *X) {

    int start_ptr = reach(L, b, buffer);

    for (int i = 0; i < b.col_begins[b.N]; i++) {
        X[b.row_indices[i]] = b.values[i];
    }

    for (int i = start_ptr; i < L.N; i++) {
        int j = buffer[i];
        if (j < 0)
            continue;

        X[j] /= L.values[L.col_begins[j]];

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
int dfs(int start_node, int start_ptr, int *buffer_upper, 
        int *buffer_lower, matrix_t<double> &L) {

    int buf_ptr = 0;
    bool is_finished = false;
    buffer_upper[buf_ptr] = start_node;
    
    while (buf_ptr >= 0) {
        is_finished = true;
        start_node = buffer_upper[buf_ptr]; 

        if (!IS_FLAGED(L.col_begins, start_node)) {
            FLAGE(L.col_begins, start_node);
            buffer_lower[buf_ptr] = (start_node < 0)? \
                    0 : RESTORE_SIGN(L.col_begins[start_node]);
        }

        int shift = (start_node < 0)? 0 : RESTORE_SIGN(L.col_begins[start_node+1]);
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
            buffer_upper[--start_ptr] = start_node;
        }
    }
    return start_ptr;
}


/**
* @brief Compute the nonzero pattern of the matrix.
* @param L: The sparse triangular matrix stored in the CSC format.
* @param b: The sparse vector stored in the CSC format.
* @param buffer: The buffer space in which the DFS operates.
* @return The beginning of the buffer space.
*/
int reach(matrix_t<double> &L, matrix_t<double> &b, int *buffer) {
    int start_ptr = L.N;
    int *buffer_upper = buffer + L.N;

    for (int i = 0; i < b.col_begins[b.N]; i++) {
        if (!IS_FLAGED(L.col_begins, b.row_indices[i])) {
            start_ptr = dfs(b.row_indices[i], start_ptr, buffer, buffer_upper, L);
        }
    }

    for (int i = start_ptr; i < L.N; i++) {
        FLAGE(L.col_begins, buffer[i]);
    }

    return start_ptr;
}


/**
* @brief The **guarded** algorithm for solving the triangular system Lx=b.
* @param dim: The matrix dimension.
* @param col_begins: The column pointer of L.
* @param row_indices: The row index of L.
* @param values: The values of L.
* @param bx: (I/O parameter) The right hand-side b at start_node 
*           and the solution x at the end.
* @return The number of FLOPS.
*/
size_t lowTriSolveGuarded(int dim, int* col_begins,
                     int* row_indices, double* values, double* bx) {
    if (!col_begins || !row_indices || !bx)
        return EMPTY_INPUT;

    size_t flop_count = 0;
    for (int j = 0; j < dim; j++) {
        if (bx[j] == 0)
            continue;

        bx[j] /= values[col_begins[j]];
        flop_count++;
        for (int p = col_begins[j] + 1; p < col_begins[j + 1]; p++) {
            bx[row_indices[p]] -= values[p] * bx[j];
            flop_count += 2;
        }
    }
    return flop_count;
}


/**
* @brief The **naive** algorithm for solving the triangular system Lx=b.
* @param dim: The matrix dimension.
* @param col_begins: The column pointer of L.
* @param row_indices: The row index of L.
* @param values: The values of L.
* @param bx: (I/O parameter) The right hand-side b at start_node 
*           and the solution x at the end.
* @return The number of FLOPS.
*/
size_t lowTriSolveNaive(int dim, int* col_begins,
                        int* row_indices, double* values, double* bx) {
    if (!col_begins || !row_indices || !bx)
        return -1;

    size_t flop_count = 0;
    for (int j = 0; j < dim; j++) {
        bx[j] /= values[col_begins[j]];
        flop_count++;
        for (int p = col_begins[j] + 1; p < col_begins[j + 1]; p++) {
            bx[row_indices[p]] -= values[p] * bx[j];
            flop_count += 2;
        }
    }
    return flop_count;
}