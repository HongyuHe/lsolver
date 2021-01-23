#include <iostream>

#include "omp.h"
#include "util.h"
#include "serial.h"
#include "verify.h"
#include "wrapper.h"


bool isResultCorrect(const matrix_t<double> &L,
                     const matrix_t<double> &b,
                     const double *const X);


int main(int argc, char **argv) {
    std::cout << "\n\tlsolver starts ...\n";
    /* File I/O */
    std::vector<std::string> fnames = getMtxFileNames(argc, argv);
    produceLowerTriMtx(fnames);

    General prop;
    matrix_t<double> L, b;
    SeldonAdapter<double> csc_L, csc_b;
    Seldon::Matrix<double, Seldon::General, 
                    Seldon::ArrayColSparse> mm_L, mm_b;

    /* CSC formatting */
    ReadMatrixMarket(fnames[1], mm_L);
    ReadMatrixMarket(fnames[2], mm_b);

    convertToCsc(mm_L, prop, csc_L);
    convertToCsc(mm_b, prop, csc_b);

    mapSparseMtx(mm_L, csc_L, L);
    mapSparseMtx(mm_b, csc_b, b);

    /* Solve the lower triangular system*/
    int X[L.N*2] = {0};
    double results[L.N] = {0};
    double wtime = omp_get_wtime();

    lowTriSolveGraph(L, b, X, results);
    writeResults(getOutputName(fnames[0]), results, L.N);

    TIMING("lowTriSolveGraph", wtime);

#ifdef TEST
    isResultCorrect(L, b, results);
#endif // TEST

    EOL;
    return EXIT_SUCCESS;
}




bool isResultCorrect(const matrix_t<double> &L,
                     const matrix_t<double> &b,
                     const double *const X) {

    double check_sum[L.M] = {0};
    for (int col = 0; col < L.N; col++) { 
        // Traverse each column.
        for (int val_ind = L.col_begins[col];
                 val_ind < L.col_begins[col+1]; 
                 val_ind++) {
            // Summing values over rows.
            int row = L.row_indices[val_ind]; 
            check_sum[row] += L.values[val_ind] * X[col];
        } 
    }
    for (int i = 0, j = 0; i < L.M && j < b.nnz; i++) {
        if (!(IS_APPROX_SAME(check_sum[i], 0.0, EPS))) { 
            // Check non-zero entries.
            if (IS_APPROX_SAME(check_sum[i], b.values[j], TOL)) { 
                debug("\t[VERIFY] Correct at b[%d]\t%.5f==%.5f\n", 
                      i, check_sum[i], b.values[j]);
            } else {
                debug("\t[VERIFY] Wrong at b[%d]\t%.5f<>%.5f\n", 
                      i, check_sum[i], b.values[j]);
                 return false;
            }
            j++;
        }
    } EOL;

    return true;
}