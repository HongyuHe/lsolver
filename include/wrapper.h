/*************************************/
/* Wrappers for external libraries   */
/*************************************/

#ifndef _WRAPPER_H
#define _WRAPPER_H

#include <iostream>
#include "SeldonFlag.hxx"
#include "SeldonLib.hxx"

#include "serial.h"

using namespace Seldon;

template <class T>
struct SeldonAdapter {
    Vector<T> values;
    IVect col_begins;
    IVect row_indices;
};

template <class T>
void mapSparseMtx (const Matrix<T, General, ArrayColSparse> &mtx,
                      const SeldonAdapter<T> &sd, matrix_t<T> &m) {
    m.M = mtx.GetM();
    m.N = mtx.GetN();
    m.nnz = sd.values.GetLength();
    m.col_begins = sd.col_begins.GetData();
    m.row_indices = sd.row_indices.GetData();
    m.values = (double*)sd.values.GetData();
}

template <class T, class Prop, class Allocator>
void printMatrix(Matrix<T, Prop, ArrayColSparse, Allocator> m) {
    cout << "[Matrix] " << m.GetM() << " X " << m.GetN() << endl; 
    for (int i = 0; i < m.GetM(); i++) {
        for (int j = 0; j < m.GetN(); j++) {
            cout << setw(5) << setprecision(3) << m.Get(i, j) << setw(5);
        } EOL;
    } EOL;
}

template <class T>
void printCsc(const struct SeldonAdapter<T>& csc) {
    cout << "[CSC format]\nValues:\n\t"; csc.values.Print();
    cout << "Row indices:\n\t"; csc.row_indices.Print();
    cout << "Col begins: \n\t"; csc.col_begins.Print();
    EOL;
}


template <class T, class Prop, class Allocator>
void getLowerTriangle(Matrix<T, Prop, ArrayColSparse, Allocator>& A,
                      Matrix<T, Prop, ArrayColSparse, Allocator>& L) {
    Copy(A, L);
    for (int i = 0; i < L.GetM(); i++) {
        for (int j = i+1; j < L.GetN(); j++) {
            if (L.Get(i, j) != 0)
                L.Set(i, j, __DBL_MIN__);
        }
    }
    L.RemoveSmallEntry(__DBL_MIN__+__DBL_EPSILON__);
}

template <class T, class Prop, class Alloc>
void convertToCsc(const Matrix<T, Prop, ArrayColSparse, Alloc> &A,
                  General &prop, struct SeldonAdapter<T>& csc) {
    ConvertToCSC(A, prop, csc.col_begins, csc.row_indices, csc.values);
    csc.col_begins.Resize(csc.col_begins.GetLength());
}

template <class T>
void clearCscBuffer(SeldonAdapter<T>& csc) {
    csc.values.Clear();
    csc.col_begins.Clear(); 
    csc.row_indices.Clear(); 
}

#endif // _WRAPPER_H