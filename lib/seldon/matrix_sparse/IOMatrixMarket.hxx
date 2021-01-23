// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2001-2011 Marc Durufle
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#ifndef SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_HXX


namespace Seldon
{
  
  template<class T>
  void ReadComplexValue(istream& FileStream, T& entry);
  
  template<class T>
  void ReadComplexValue(istream& FileStream, complex<T>& entry);

  template<class T>
  void WriteComplexValue(ostream& FileStream, const T& entry);

  template<class T>
  void WriteComplexValue(ostream& FileStream, const complex<T>& entry);

  template<class Tint, class AllocInt, class T, class Allocator>
  void ReadCoordinateMatrix(istream& FileStream,
                            Vector<Tint, VectFull, AllocInt>& row_numbers,
                            Vector<Tint, VectFull, AllocInt>& col_numbers,
                            Vector<T, VectFull, Allocator>& values,
                            bool cplx = false);
  
  template<class Matrix1, class T>
  void ReadCoordinateMatrix(Matrix1& A, istream& FileStream, T& zero,
                            int index = 1, long nnz = -1, bool cplx = false);

  template<class Tint, class AllocInt, class T, class Allocator>
  void WriteCoordinateMatrix(ostream& FileStream,
			     const Vector<Tint, VectFull, AllocInt>& row_numbers,
			     const Vector<Tint, VectFull, AllocInt>& col_numbers,
			     const Vector<T, VectFull, Allocator>& values,
			     bool cplx = false);
  
  template<class T0, class Prop0, class Storage0, class Alloc0, class T>
  void WriteCoordinateMatrix(const Matrix<T0, Prop0, Storage0, Alloc0>& A,
                             ostream& FileStream, T& zero,
                             int index = 1, bool cplx = false);

  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<T, Prop, ColSparse, Allocator>& A);


  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<complex<T>, Prop, ColSparse, Allocator>& A);


  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<T, Prop, RowSymSparse, Allocator>& A);


  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<complex<T>, Prop, RowSymSparse, Allocator>& A);


  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const T2& val, 
                         Matrix<T, Symmetric, Storage, Allocator>& A);


  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const T2& val, 
                         Matrix<T, General, Storage, Allocator>& A);


  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const complex<T2>& val, 
                         Matrix<T, Symmetric, Storage, Allocator>& A);


  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const complex<T2>& val, 
                         Matrix<T, General, Storage, Allocator>& A);


  template <class T, class Prop, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<T, Prop, Storage, Allocator>& A);


  template<class T>
  void ReadComplexValuesHarwell(long Nnonzero, int Nline_val, int line_width_val,
                                int element_width_val,
                                istream& input_stream, T* A_data);
  
  
  template<class T>
  void ReadComplexValuesHarwell(long Nnonzero, int Nline_val, int line_width_val,
                                int element_width_val,
                                istream& input_stream, complex<T>* A_data);
  
  
  template <class T, class Prop, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename,
                         string value_type, string matrix_type,
                         Matrix<T, Prop, Storage, Allocator> & A);
  
  
  template<class T>
  void PrintComplexValuesHarwell(long nnz, const Vector<complex<T> >& Val,
                                 ofstream& file_out);


  template<class T>
  void PrintComplexValuesHarwell(long nnz, const Vector<T>& Val,
                                 ofstream& file_out);


  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<T, Prop, Storage, Allocator>& A,
                          const string& file_name);
  

  template <class T, class Prop, class Storage, class Allocator>
  void ReadMatrixMarket(string filename,
                        Matrix<T, Prop, Storage, Allocator>& A);
  

  template<class T, class Prop, class Storage, class Allocator>
  void WriteMatrixMarket(const Matrix<T, Prop, Storage, Allocator>& A,
                         const string& file_name);
  
  
}


#define SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_HXX
#endif
