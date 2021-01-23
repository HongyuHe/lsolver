#ifndef _UTIL_H
#define _UTIL_H

#include <iostream>
#include <iomanip>
#include <string>
#include <limits>
#include <vector>
#include <cmath>
#include <omp.h>

#define EMPTY_INPUT -1
#define EPS pow(10, -12)
#define TOL pow(10, -10)

#define IS_APPROX_SAME(a, b, tol) abs(a - b) < tol

#define TIMING(x, t)  do { \
                std::cout << std::boolalpha \
                << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "\t[TIMING] " \
                << x << ":\t" << std::fixed << (omp_get_wtime() - t)*1000 << " ms"<< std::endl; } while (0)

#define DEBUG(x)  do { \
                std::cerr << std::boolalpha \
                << std::setprecision(std::numeric_limits<double>::digits10 + 1) << "\t[DEBUG] " \
                << x << std::endl; } while (0)
#define DEBUG2(x) do { std::cerr << "\t[DEBUG] "<< #x << ": " << x << std::endl; } while (0)


#define debug(...) \
            fprintf(stderr,  __VA_ARGS__)
            
#define EOL std::cout << std::endl


std::vector<std::string> getMtxFileNames(int argc, char **argv);

void produceLowerTriMtx(std::vector<std::string> fns);

void printVector(const double *const vec, int n);

std::string getOutputName(std::string fname);

bool fileExists (const std::string& name);

void writeResults(const std::string& fname, 
                 const double *const X, int dim);
                 
#endif // _UTIL_H

