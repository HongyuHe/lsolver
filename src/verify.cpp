#include <iostream>

#include "serial.h"
#include "util.h"

bool isResultCorrect(const matrix_t<double> *const L,
                     const matrix_t<double> *const b,
                     const double *const X) {
    double check_sum[L->M] = {0};
    for (int col = 0; col < L->N; col++) { 
        // Traverse each column.
        for (int val_ind = L->col_begins[col];
                 val_ind < L->col_begins[col+1]; 
                 val_ind++) {             
            // Summing values over rows.
            int row = L->row_indices[val_ind];
            check_sum[row] += L->values[val_ind] * X[col];
        } 
    }
    for (int i = 0, j = 0; i < L->M && j < b->nnz; i++) {
        if (!(IS_APPROX_SAME(check_sum[i], 0.0, EPS))) { 
            // Check non-zero entries.
            if ((IS_APPROX_SAME(check_sum[i], b->values[j], TOL))) {
                return false;
            }
            j++;
        }
    }
    return true;
}