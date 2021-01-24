#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vector>
#include "mm_io.h"

#include "util.h"

std::vector<std::string> getMtxFileNames(int argc, char **argv) {
    std::string fn;
    std::vector<std::string> fnames;

    if (argc < 3) {
        fprintf(stderr, 
                "Usage: %s [martix-market-file-of-L] [martix-market-file-of-b]\n", argv[0]);
        exit(EXIT_FAILURE);
    } else {
        fn = argv[1];
        fnames.push_back(argv[1]);
        fnames.push_back(fn.substr(0, fn.find(".mtx"))+"_tri.mtx");
        fnames.push_back(argv[2]);
    }
    return fnames;
}

void produceLowerTriMtx(std::vector<std::string> fnames) {
    MM_typecode matcode;
    FILE *fh_in, *fh_out;
    int M, N, count;
    int i, *I, *J;
    int ret_code;
    double *val;

    if ((fh_in = fopen(fnames[0].c_str(), "r")) == NULL || 
       ((fh_out = fopen(fnames[1].c_str(), "w")) == NULL)) {
           DEBUG("IO failed");
           exit(EXIT_FAILURE);
       }
    if (mm_read_banner(fh_in, &matcode) != 0) {
        DEBUG("No Matrix Market banner.\n");
        exit(EXIT_FAILURE);
    }
    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode)) {
        DEBUG("No support.\n");
        exit(EXIT_FAILURE);
    }
    if ((ret_code = mm_read_mtx_crd_size(fh_in, &M, &N, &count)) != 0)
        exit(EXIT_FAILURE);

    /* Reseve memory for matrices */

    I = (int *)malloc(count * sizeof(int));
    J = (int *)malloc(count * sizeof(int));
    val = (double *)malloc(count * sizeof(double));

    int tri_count = 0;
    for (i = 0; i < count; i++) {
        fscanf(fh_in, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        if (I[i] >= J[i]) 
            tri_count++;
    }

    if (fh_in != stdin)
        fclose(fh_in);

    mm_write_banner(fh_out, matcode);
    mm_write_mtx_crd_size(fh_out, M, N, tri_count);

    for (i = 0; i < count; i++) {
        if (I[i] >= J[i]) {
            fprintf(fh_out, "%d %d %20.16g\n", I[i], J[i], val[i]);
        }
    }

    if (fh_in != stdin)
        fclose(fh_out);
}

std::string getOutputName(std::string fname) {
    return fname.substr(0, fname.find(".mtx"))
                .replace(fname.find("matrices"), 
                         sizeof("matrices") - 1, "out");
}

void writeResults(const std::string& fname, 
                 const double *const X, int dim) {
    std::cout << "\n\t[OUTPUT] Output file:\t\t" << fname;
    std::ofstream fh;
    fh.open(fname);

    for (int i = 0; i < dim; i++) {
        fh << X[i] << std::endl;
    }
    fh.close();
    EOL;
}

bool fileExists (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

void printVector(const double *const vec, int n) {
    for (int i = 0; i < n; i++) {
        std::cout << vec[i] << "\t";
    } EOL;
}
