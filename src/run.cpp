#include <iostream>

#include "omp.h"
#include "util.h"
#include "serial.h"
#include "verify.h"
#include "wrapper.h"


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
    isResultCorrect(&L, &b, results);
#endif // TEST

    EOL;
    return EXIT_SUCCESS;
}