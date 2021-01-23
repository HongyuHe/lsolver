#include <iostream>

#include "util.h"
#include "serial.h"
#include "verify.h"
#include "wrapper.h"

int main() {
    std::cout << "\n\t[TEST 2] trivial2.mtx\n";

    General prop;
    SeldonAdapter<double> csc_L, csc_b;
    Seldon::Matrix<double, Seldon::General, Seldon::ArrayColSparse> mm_L, mm_b;
    std::vector<std::string> fnames = {"./data/matrices/trivial/trivial2.mtx",
                                      "./data/matrices/trivial/trivial2_tri.mtx",
                                      "./data/matrices/trivial/trivial2_b.mtx"};

    double wtime = omp_get_wtime();
    produceLowerTriMtx(fnames);
    ReadMatrixMarket(fnames[1], mm_L);
    ReadMatrixMarket(fnames[2], mm_b);
    convertToCsc(mm_L, prop, csc_L);
    convertToCsc(mm_b, prop, csc_b);
    TIMING("I/O & CSC conversion", wtime);
    EOL;

    matrix_t<double> L, b;
    mapSparseMtx(mm_L, csc_L, L);
    mapSparseMtx(mm_b, csc_b, b);

    double bx[L.N] = {0};
    for (int i = 0; i < b.nnz; i++) {
        bx[b.row_indices[i]] = b.values[i]; 
    }
    wtime = omp_get_wtime();
    lowTriSolveNaive(L.N, L.col_begins, L.row_indices, L.values, bx);
    TIMING("lowTriSolveNaive", wtime);

    for (int i = 0; i < b.nnz; i++) {
        bx[b.row_indices[i]] = b.values[i]; 
    }
    wtime = omp_get_wtime();
    lowTriSolveGuarded(L.N, L.col_begins, L.row_indices, L.values, bx);
    TIMING("lowTriSolveGuarded", wtime);

    int X[L.N*2] = {0};
    double results[L.N] = {0};
    wtime = omp_get_wtime();
    lowTriSolveGraph(L, b, X, results);
    TIMING("lowTriSolveGraph", wtime);

    EOL;
    return EXIT_SUCCESS;
}