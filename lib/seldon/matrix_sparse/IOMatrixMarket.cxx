// Copyright (C) 2001-2011 Vivien Mallet
// Copyright (C) 2001-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_CXX


#include "IOMatrixMarket.hxx"
#include <iomanip>

/*
  Functions defined in this file:
  (storage RowSparse, ColSparse, RowSymSparse, ColSymSparse,
  ArrayRowSparse, ArrayColSparse, ArrayRowSymSparse, ArrayColSymSparse
  and equivalent complex storages if SeldonComplexMatrix.hxx is included)
  
  A is read in Harwell-Boeing format (.rua, .rsa, .rra, .cua, .cra or .csa)
  ReadHarwellBoeing(file_name, A)
  
  A is written in Harwell-Boeing format
  WriteHarwellBoeing(A, file_name)
  
  A is read in Matrix-Market format (.mtx)
  ReadMatrixMarket(file_name, A)

  A is written in Matrix-Market format (.mtx)
  WriteMatrixMarket(A, file_name)

  Internal functions (do not use directly) :
  ReadComplexValuesHarwell
  PrintComplexValuesHarwell
  
*/

namespace Seldon
{


  //////////////////////////
  // ReadCoordinateMatrix //
  //////////////////////////
  

  //! reads a real or complex value in a file
  template<class T>
  void ReadComplexValue(istream& FileStream, T& entry)
  {
    FileStream >> entry;
  }
  
  
  //! reads a real or complex value in a file
  template<class T>
  void ReadComplexValue(istream& FileStream, complex<T>& entry)
  {
    T a, b;
    FileStream >> a >> b;
    entry = complex<T>(a, b);
  }
  

  //! reads a real or complex value in a file
  template<class T>
  void WriteComplexValue(ostream& FileStream, const T& entry)
  {
    FileStream << entry;
  }
  
  
  //! reads a real or complex value in a file
  template<class T>
  void WriteComplexValue(ostream& FileStream, const complex<T>& entry)
  {
    FileStream << real(entry) << " " << imag(entry);
  }
  
    
  //! Reading of matrix in coordinate format
  /*!
    \param FileStream stream where the matrix is read
    \param row_numbers row indices
    \param col_numbers column indices
    \param values
    \param cplx if true, imaginary and real part are in
                two separate columns
   */  
  template<class Tint, class AllocInt, class T, class Allocator>
  void ReadCoordinateMatrix(istream& FileStream,
                            Vector<Tint, VectFull, AllocInt>& row_numbers,
                            Vector<Tint, VectFull, AllocInt>& col_numbers,
                            Vector<T, VectFull, Allocator>& values,
                            bool cplx)
  {
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("ReadCoordinateMatrix", "Stream is not ready.");
#endif

    T entry; int row = 0, col = 0;
    long nb_elt = 0;
    SetComplexZero(entry);

    while (!FileStream.eof())
      {
	// new entry is read (1-index)
	FileStream >> row >> col;
        if (cplx)
          ReadComplexValue(FileStream, entry);
        else
          FileStream >> entry;
        
	if (FileStream.fail())
	  break;
	else
	  {
#ifdef SELDON_CHECK_IO
	    if (row < 1)
	      throw IOError(string("ReadCoordinateMatrix"),
			    string("Error : Row number should be greater ")
			    + "than 0 but is equal to " + to_str(row));

	    if (col < 1)
	      throw IOError(string("ReadCoordinateMatrix"),
			    string("Error : Column number should be greater")
			    + " than 0 but is equal to " + to_str(col));
#endif

	    nb_elt++;

	    // inserting new element
	    if (nb_elt > values.GetM())
	      {
		values.Resize(2*nb_elt);
		row_numbers.Resize(2*nb_elt);
		col_numbers.Resize(2*nb_elt);
	      }

	    values(nb_elt-1) = entry;
	    row_numbers(nb_elt-1) = row;
	    col_numbers(nb_elt-1) = col;
	  }
      }
    
    if (nb_elt > 0)
      {
	row_numbers.Resize(nb_elt);
	col_numbers.Resize(nb_elt);
	values.Resize(nb_elt);
      }
  }
  

  //! Reading of matrix in coordinate format
  /*!
    \param A output matrix
    \param FileStream stream where the matrix is read
    \param zero type of value to read (double or complex<double>)
    \param index starting index (usually 0 or 1)
    \param nnz number of non-zero entries
    \param cplx if true, imaginary and real part are in
                two separate columns
    If nnz is equal to -1, we consider that the number of non-zero entries
    is unknown and is deduced from the number of lines present in the stream
   */
  template<class Matrix1, class T>
  void ReadCoordinateMatrix(Matrix1& A, istream& FileStream, T& zero,
                            int index, long nnz, bool cplx)
  {
    // previous elements are removed
    A.Clear();
    
    Vector<int> row_numbers, col_numbers;
    Vector<T> values;
    if (nnz >= 0)
      {
        values.Reallocate(nnz);
        row_numbers.Reallocate(nnz);
        col_numbers.Reallocate(nnz);
	values.Fill(0);
	row_numbers.Zero();
	col_numbers.Zero();
      }
    
    ReadCoordinateMatrix(FileStream, row_numbers, col_numbers, values, cplx);
    
    if (row_numbers.GetM() > 0)
      ConvertMatrix_from_Coordinates(row_numbers, col_numbers, values,
                                     A, index);
  }


  //! Writes a matrix in coordinate format
  /*!
    \param FileStream stream where the matrix is read
    \param row_numbers row indices
    \param col_numbers column indices
    \param values
    \param cplx if true, imaginary and real part are written in
                two separate columns
   */    
  template<class Tint, class AllocInt, class T, class Allocator>
  void WriteCoordinateMatrix(ostream& FileStream,
                             const Vector<Tint, VectFull, AllocInt>& row_numbers,
                             const Vector<Tint, VectFull, AllocInt>& col_numbers,
                             const Vector<T, VectFull, Allocator>& values,
                             bool cplx)
  {
    if (cplx)
      {
        for (int i = 0; i < row_numbers.GetM(); i++)
          {
            FileStream << row_numbers(i) << " " << col_numbers(i) << " ";
            WriteComplexValue(FileStream, values(i));
            FileStream << '\n';    
          }
      }
    else
      {
        for (int i = 0; i < row_numbers.GetM(); i++)
          FileStream << row_numbers(i) << " " << col_numbers(i)
                     << " " << values(i) << '\n';
      }
  }
  
  
  //! Writes matrix in coordinate format
  /*!
    \param A output matrix
    \param FileStream stream where the matrix is read
    \param zero type of value to read (double or complex<double>)
    \param index starting index (usually 0 or 1)
    \param cplx if true, imaginary and real part are in
                two separate columns
  */
  template<class T0, class Prop0, class Storage0, class Alloc0, class T>
  void WriteCoordinateMatrix(const Matrix<T0, Prop0, Storage0, Alloc0>& A,
                             ostream& FileStream, T& zero,
                             int index, bool cplx)
  {
    // conversion to coordinate format (if symmetric part, lower and upper part
    // are recovered)
    Vector<int> IndRow, IndCol;
    Vector<T> Value;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol,
                                 Value, index, true);
    
    WriteCoordinateMatrix(FileStream, IndRow, IndCol, Value, cplx);
    
    // if last element is not present, 0 is printed
    long N = IndRow.GetM()-1;
    if (N >= 0)
      {
        int m = A.GetM()-1, n = A.GetN()-1;
        SetComplexZero(zero);
        if ( (IndRow(N) != m+index) || (IndCol(N) != n+index))
          {
            if (A(m, n) == zero)
              {
                FileStream << m+index << " " << n+index << " ";
                if (cplx)
                  WriteComplexValue(FileStream, zero);
                else
                  FileStream << zero;
                
                FileStream << '\n';
              }          
          }
      }
  }
  
  
  //! Reads a sparse matrix in a file in Harwell-Boeing format.
  /*! This functions was written to read the files in Harwell-Boeing format as
    distributed on the Matrix Market, http://math.nist.gov/MatrixMarket/. A
    file from Matrix Market is associated with a type encoded in three
    characters (which is usually its extension). The supported types have:
    - 'R' (real), 'P' (pattern), or 'C' (complex) as first character
    - 'U' (unsymmetric), 'S' (symmetric) or 'R' (rectangular) as second character;
    'H' (Hermitian) and 'Z' (skew symmetric) are not supported;
    - 'A' (assembled) as third character; 'E' (elemental) is not supported.
    \param[in] filename path to the file that contains the matrix.
    \param[out] A the matrix to be structured and filled with the values of \a
    filename.
  */
  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<T, Prop, ColSparse, Allocator>& A)
  {
    ReadHarwellBoeing(filename, "real", "general", A);
  }


  //! \copydoc ReadHarwellBoeing(string filename,
  //!  Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<complex<T>, Prop, ColSparse, Allocator>& A)
  {
    ReadHarwellBoeing(filename, "complex", "general", A);
  }

  
  //! \copydoc ReadHarwellBoeing(string filename,
  //! Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<T, Prop, RowSymSparse, Allocator>& A)
  {
    ReadHarwellBoeing(filename, "real", "symmetric", A);
  }

  
  //! \copydoc ReadHarwellBoeing(string filename,
  //! Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class Prop, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<complex<T>, Prop, RowSymSparse, Allocator>& A)
  {
    ReadHarwellBoeing(filename, "complex", "symmetric", A);
  }

  
  //! \copydoc ReadHarwellBoeing(string filename,
  //! Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const T2& val, 
                         Matrix<T, Symmetric, Storage, Allocator>& A)
  {  
    Matrix<T2, Symmetric, RowSymSparse> B;
    ReadHarwellBoeing(filename, "real", "symmetric", B);       
    Copy(B, A);
  }
  
  
  //! \copydoc ReadHarwellBoeing(string filename,
  //! Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const T2& val, 
                         Matrix<T, General, Storage, Allocator>& A)
  {
    Matrix<T2, General, ColSparse> B;
    ReadHarwellBoeing(filename, "real", "general", B);
    Copy(B, A);
  }
  
  
  //! \copydoc ReadHarwellBoeing(string filename,
  //! Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const complex<T2>& val, 
                         Matrix<T, Symmetric, Storage, Allocator>& A)
  {  
    Matrix<complex<T2>, Symmetric, RowSymSparse> B;
    ReadHarwellBoeing(filename, "complex", "symmetric", B);       
    Copy(B, A);
  }
  
  
  //! \copydoc ReadHarwellBoeing(string filename,
  //! Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class T2, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename, const complex<T2>& val, 
                         Matrix<T, General, Storage, Allocator>& A)
  {
    Matrix<complex<T2>, General, ColSparse> B;
    ReadHarwellBoeing(filename, "complex", "general", B);
    Copy(B, A);
  }
  
  
  //! \copydoc ReadHarwellBoeing(string filename,
  //!  Matrix<T, Prop, ColSparse, Allocator>& A)
  template <class T, class Prop, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename,
                         Matrix<T, Prop, Storage, Allocator>& A)
  {
    typename Matrix<T, Prop, Storage, Allocator>::entry_type val;
    ReadHarwellBoeing(filename, val, A);
  }
  
  
  template<class T>
  void ReadComplexValuesHarwell(long Nnonzero, int Nline_val, int line_width_val,
                                int element_width_val,
                                istream& input_stream, T* A_data)
  {
    string line, value;
    string::size_type pos;
    long index = 0;
    for (int i = 0; i < Nline_val; i++)
      {
        getline(input_stream, line);
        int k = 0;
        for (int j = 0; j < line_width_val; j++)
          {
            value = line.substr(k, element_width_val);
            pos = value.find('D');
            if (pos != string::npos)
              value[pos] = 'E';
            
            istringstream flux(value);
            flux >> A_data[index];
            index++;
            if (index == Nnonzero)
              // So as not to read more elements than actually available
              // on the line.
              break;
            k += element_width_val;
          }
        if (index == Nnonzero)
          break;
      }
  }
  
  
  template<class T>
  void ReadComplexValuesHarwell(long Nnonzero, int Nline_val, int line_width_val,
                                int element_width_val,
                                istream& input_stream, complex<T>* A_data)
  {
    string line, value;
    string::size_type pos;
    long index = 0, nb = 0; T a, b(0);
    for (int i = 0; i < Nline_val; i++)
      {
        getline(input_stream, line);
        int k = 0;
        for (int j = 0; j < line_width_val; j++)
          {
            a = b;
            value = line.substr(k, element_width_val);
            pos = value.find('D');
            if (pos != string::npos)
              value[pos] = 'E';
            
            istringstream flux(value);
            flux >> b;
            if (index%2 == 1)
              {
                A_data[nb] = complex<T>(a, b);
                nb++;
              }
            
            index++;
            if (index == 2*Nnonzero)
              // So as not to read more elements than actually available
              // on the line.
              break;
            k += element_width_val;
          }
        if (index == 2*Nnonzero)
          break;
      }
  }
  
  
  //! Reads a sparse matrix in a file in Harwell-Boeing format.
  /*! This functions was written to read the files in Harwell-Boeing format as
    distributed on the Matrix Market, http://math.nist.gov/MatrixMarket/. A
    file from Matrix Market is associated with a type encoded in three
    characters (which is usually its extension). The supported types have:
    - 'R' (real), 'P' (pattern) or 'C' (complex) as first character;
    - 'U' (unsymmetric), 'S' (symmetric) or 'R' (rectangular) as second character;
     'H' (Hermitian) and 'Z' (skew symmetric) are not supported;
    - 'A' (assembled) as third character; 'E' (elemental) is not supported.
    \param[in] filename path to the file that contains the matrix.
    \param[in] value_type the type of involved data: only "real" is supported.
    \param[in] matrix_type the type of matrix: only "general" is supported.
    \param[out] A the matrix to be structured and filled with the values of \a
    filename.
    \note This function is not supposed to be called directly unless the user
    knows what he/she is doing. See overloaded functions ReadHarwellBoeing
    with only \a filename and \a A as arguments.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void ReadHarwellBoeing(string filename,
                         string value_type, string matrix_type,
                         Matrix<T, Prop, Storage, Allocator> & A)
  {
    int i, j, k;
    string line, element;

    // Dimensions of the output matrix.
    int Nrow, Ncol;
    long Nnonzero;

    ifstream input_stream(filename.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!input_stream.good())
      throw IOError("ReadHarwellBoeing(string filename, "
                    "Matrix& A, string value_type, string matrix_type)",
                    "Unable to read file \"" + filename + "\".");
#endif

    /*** Parsing headers ***/

    // First line.
    string title, name;
    getline(input_stream, line);
    title = line.substr(0, 72);
    name = line.substr(72, 8);

    // Second line.
    int Nline_ptr, Nline_ind, Nline_val, Nline_rhs;
    input_stream >> i >> Nline_ptr >> Nline_ind >> Nline_val >> Nline_rhs;
    
    // Third line.
    string type;
    int Nelemental;
    input_stream >> type >> Nrow >> Ncol >> Nnonzero >> Nelemental;

#ifdef SELDON_CHECK_IO
    // Checks if the stream is still in good shape.
    if (!input_stream.good())
      throw IOError("ReadHarwellBoeing(string filename, "
                    "Matrix& A, string value_type, string matrix_type)",
                    "Unable to read integers in file \""
                    + filename + "\".");
#endif

    if (type.size() != 3)
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a matrix of "
                          "type \"" + type + "\", which is not a valid type. "
                          "The type must contain exactly three characters.");

    if (type.substr(0, 1) != "R" && type.substr(0, 1) != "C"
        && type.substr(0, 1) != "P")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a matrix of "
                          "type \"" + type + "\", which is not a valid type. "
                          "The first character of that type must be 'R', 'C' "
                          "(not supported anyway) or 'P'.");

    if (type.substr(0, 1) == "R" && value_type != "real")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a real-valued "
                          "matrix, while the input matrix 'A' is not "
                          "declared as such.");

    if (type.substr(0, 1) == "C" && value_type != "complex")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a complex-valued "
                          "matrix, while the input matrix 'A' is not "
                          "declared as such.");
    
    if (type.substr(1, 1) == "H")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a Hermitian "
                          "matrix, which is not supported.");

    if (type.substr(1, 1) == "Z")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains a skew "
                          "symmetric matrix, which is not supported.");
    
    if (type.substr(2, 1) == "E")
      throw WrongArgument("ReadHarwellBoeing(string filename, "
                          "Matrix& A, string value_type, string matrix_type)",
                          "File \"" + filename + "\" contains an elemental "
                          "matrix, which is not supported.");

    getline(input_stream, line);

    // Fourth line.
    getline(input_stream, line);
    int line_width_ptr(0), element_width_ptr(0);
    int line_width_ind(0), element_width_ind(0);
    string::size_type pos1 = line.find('(', 0);
    string::size_type pos2 = line.find('I', 0);
    if ((pos2 < 16) && (pos1 < pos2))
      line_width_ptr = to_num<int>(line.substr(pos1+1, pos2-pos1-1));
    
    pos1 = pos2;
    pos2 = line.find(')', 0);
    if ((pos2 < 16) && (pos1 < pos2))
      element_width_ptr = to_num<int>(line.substr(pos1+1, pos2-pos1-1));

    pos1 = line.find('(', 16);
    pos2 = line.find('I', 16);
    if ((pos2 < 32) && (pos1 < pos2))
      line_width_ind = to_num<int>(line.substr(pos1+1, pos2-pos1-1));

    pos1 = pos2;
    pos2 = line.find(')', 16);
    if ((pos2 < 32) && (pos1 < pos2))
      element_width_ind = to_num<int>(line.substr(pos1+1, pos2-pos1-1));
    
    int line_width_val = 0;
    int element_width_val = 0;
    if (type.substr(0, 1) != "P")
      {
        element = line.substr(32, 20);

        // Splits the format to retrieve the useful widths. This part is
        // tricky because no parsing (in C++) of Fortran format was found.
        Vector<int> vect;
        string delimiter = " (ABDEFGILOPZ*.)";
        string tmp_str;
        int tmp_int;
        string::size_type index_beg = element.find_first_not_of(delimiter);
        string::size_type index_end;
        while (index_beg != string::npos)
          {
            index_end = element.find_first_of(delimiter, index_beg);
            tmp_str = element.substr(index_beg,
                                     index_end == string::npos ?
                                     string::npos : (index_end-index_beg));
            istringstream(tmp_str) >> tmp_int;
            vect.PushBack(tmp_int);
            index_beg = element.find_first_not_of(delimiter, index_end);
          }
        
        if (vect.GetM() < 3)
          throw WrongArgument("ReadHarwellBoeing(string filename, Matrix& A, "
                              "string value_type, string matrix_type)",
                              "File \"" + filename + "\" contains values "
                              "written in format \"" + element + "\", which "
                              "could not be parsed.");

        line_width_val = vect(vect.GetM() - 3);
        element_width_val = vect(vect.GetM() - 2);
      }

    if ((line_width_ptr == 0) || (element_width_ptr == 0) ||
        (line_width_ind == 0) || (element_width_ind == 0) ||
        (line_width_val == 0) || (element_width_val == 0) )
      throw WrongArgument("ReadHarwellBoeing(string filename, Matrix& A, "
                          "string value_type, string matrix_type)",
                          "Failed to read formats in file"+filename);
    
    // Fifth header line, if any: ignored. RHS are not read.
    if (Nline_rhs != 0)
      getline(input_stream, line);

    /*** Allocations ***/
    
    typedef typename SeldonDefaultAllocator<VectFull, int>::allocator
      AllocatorInt;

    typedef typename SeldonDefaultAllocator<VectFull, long>::allocator
      AllocatorLong;
    
    // Content of output matrix A.
    long* A_ptr;
    int* A_ind;
    T* A_data;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	A_ptr = reinterpret_cast<long*>(AllocatorLong::
                                        allocate(Ncol + 1));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        A_ptr = NULL;
      }

    if (A_ptr == NULL)
      throw NoMemory("ReadHarwellBoeing(string filename, "
                     "Matrix& A, string value_type, string matrix_type)",
		     "Unable to allocate memory for an array of "
		     + to_str(Ncol + 1) + " integers.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

        // Reallocates 'A_ind' and 'A_data' in order to append the
        // elements of the i-th row of C.
        A_ind = reinterpret_cast<int*>(AllocatorInt::allocate(Nnonzero));
        A_data = reinterpret_cast<T*>
          (Allocator::allocate(Nnonzero));

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
        A_ind = NULL;
        A_data = NULL;
      }

    if (A_ind == NULL || A_data == NULL)
      throw NoMemory("ReadHarwellBoeing(string filename, "
                     "Matrix& A, string value_type, string matrix_type)",
                     "Unable to allocate memory for an array of "
                     + to_str(Nnonzero) + " integers "
                     "and for an array of "
                     + to_str(sizeof(T) * Nnonzero) + " bytes.");
#endif

    /*** Reads the structure ***/

    long index = 0;
    for (i = 0; i < Nline_ptr; i++)
      {
        getline(input_stream, line);
        k = 0;
        for (j = 0; j < line_width_ptr; j++)
          {
            istringstream(line.substr(k, element_width_ptr)) >> A_ptr[index];
            
            // The indexes are 1-based, so this corrects it:
            A_ptr[index]--;
            index++;
            if (index == Ncol + 1)
              // So as not to read more elements than actually available on
              // the line.
              break;
            k += element_width_ptr;
          }
        if (index == Ncol + 1)
          break;
      }

    index = 0;
    for (i = 0; i < Nline_ind; i++)
      {
        getline(input_stream, line);
        k = 0;
        for (j = 0; j < line_width_ind; j++)
          {
            istringstream(line.substr(k, element_width_ind)) >> A_ind[index];
            // The indexes are 1-based, so this corrects it:
            A_ind[index]--;
            index++;
            if (index == Nnonzero)
              // So as not to read more elements than actually available on
              // the line.
              break;
            k += element_width_ind;
          }
        if (index == Nnonzero)
          break;
      }
    
    /*** Reads the values ***/

    if (type.substr(0, 1) == "P")
      for (long i0 = 0; i0 < Nnonzero; i0++)
        A_data[i0] = T(1);
    else
      {
        index = 0;
        ReadComplexValuesHarwell(Nnonzero, Nline_val, line_width_val,
                                 element_width_val, input_stream, A_data);
      }

#ifdef SELDON_CHECK_IO
    // Checks if the stream is still in good shape.
    if (!input_stream.good())
      throw IOError("ReadHarwellBoeing(string filename, "
                    "Matrix& A, string value_type, string matrix_type)",
                    "Unable to read all values in file \""
                    + filename + "\".");
#endif

    input_stream.close();

    A.SetData(Nrow, Ncol, Nnonzero, A_data, A_ptr, A_ind);
  }

  
  template<class T>
  void PrintComplexValuesHarwell(long nnz, const Vector<complex<T> >& Val,
                                 ofstream& file_out)
  {
    for (long i = 0; i < 2*nnz; i += 3)
      {
        for (long j = i; j < min(i+3, 2*nnz); j++)
          {
            if (j%2 == 0)
              file_out << setw(23) << std::real(Val(j/2));
            else
              file_out << setw(23) << std::imag(Val(j/2));
          }
        
        file_out << '\n';
      }
  }


  template<class T>
  void PrintComplexValuesHarwell(long nnz, const Vector<T>& Val,
                                 ofstream& file_out)
  {
    for (long i = 0; i < nnz; i += 3)
      {
        for (long j = i; j < min(i+3, nnz); j++)
          file_out << setw(23) << Val(j);
        
        file_out << '\n';
      }
  }
  
  
  //! A is written in Harwell-Boeing file
  template<class T, class Prop, class Storage, class Allocator>
  void WriteHarwellBoeing(const Matrix<T, Prop, Storage, Allocator>& A,
                          const string& file_name)
  {
    typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
    bool complex = ( sizeof(Tcplx)/sizeof(Treal) == 2);
    
    // converting to CSC
    Vector<long> Ptr; Vector<int> Ind;
    Vector<Tcplx> Val;
    Prop sym;
    if (IsSymmetricMatrix(A))
      ConvertToCSR(A, sym, Ptr, Ind, Val);
    else
      ConvertToCSC(A, sym, Ptr, Ind, Val);
    
    // number of columns and non-zero entries
    long nnz = Val.GetM();
    int m = A.GetM();
    int n = A.GetN();
    
    // First line
    ofstream file_out(file_name.data());

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!file_out.good())
      throw IOError("WriteHarwellBoeing(const Matrix& A, string filename)",
                    "Unable to open file \"" + file_name + "\".");
#endif
    
    string title("Seldon");
    string key("S0000008");
    file_out.setf(ios::left);
    file_out << setw(72) << title;
    file_out << setw(8) << key;
    file_out << endl;

    // Compute column pointer format
    int ptr_len = int(ceil( log10(0.1 + nnz + 1) )) + 1;
    int ptr_nperline = min(80 / ptr_len, n+1);
    int ptrcrd = n / ptr_nperline + 1;
    string ptrfmt = string("(") + to_str(ptr_nperline) 
      + "I" + to_str(ptr_len) + ")";
    
    // Compute row index format
    int ind_len = int(ceil( log10(0.1 + n + 1) )) + 1;
    int ind_nperline = min(80 / ind_len, int(nnz));
    int indcrd = (nnz-1) / ind_nperline + 1;
    string indfmt = string("(") + to_str(ind_nperline)
      + 'I' + to_str(ind_len) + ')';

    // compute value format
    string valfmt("(3D23.15)");
    int valcrd = (nnz-1) / 3 + 1;
    if (complex)
      valcrd = (2*nnz-1) / 3 + 1;
    
    int rhscrd = 0;
    int totcrd = ptrcrd + indcrd + valcrd + rhscrd;
    
    // Second line
    file_out << setw(14) << totcrd;
    file_out << setw(14) << ptrcrd;
    file_out << setw(14) << indcrd;
    file_out << setw(14) << valcrd;
    file_out << setw(14) << rhscrd;
    file_out << endl;

    // Third line
    int neltvl = 0;
    char mxtype[4];
    if (complex)
      mxtype[0] = 'C';
    else
      mxtype[0] = 'R';
    
    if (m == n)
      {
        if (IsSymmetricMatrix(A))
          mxtype[1] = 'S';
        else
          mxtype[1] = 'U';
      }
    else
      mxtype[1] = 'R';
    
    mxtype[2] = 'A';
    mxtype[3] = '\0';      
    
    file_out << mxtype << "           ";
    file_out << setw(14) << m;
    file_out << setw(14) << n;
    file_out << setw(14) << nnz;
    file_out << setw(14) << neltvl;
    file_out << endl;

    // Fourth line
    file_out << setw(16) << ptrfmt;
    file_out << setw(16) << indfmt;
    file_out << setw(20) << valfmt;
    file_out << setw(20) << valfmt;
    file_out << endl;

    // printing ptr values
    file_out.setf(ios::left);
    for (int i = 0; i <= n; i += ptr_nperline)
      {
        for (int j = i; j < min(i+ptr_nperline, n+1); j++)
          file_out << setw(ptr_len) << Ptr(j)+1;

        file_out << '\n';
      }

    // printing ind values
    for (long i = 0; i < nnz; i += ind_nperline)
      {
        for (long j = i; j < min(i+ind_nperline, nnz); j++)
          file_out << setw(ind_len) << Ind(j)+1;
        
        file_out << '\n';
      }
    
    // printing values
    file_out.precision(15);
    file_out.setf(ios::scientific);
    PrintComplexValuesHarwell(nnz, Val, file_out);
    
    file_out.close();
  }

  
  //! Reads a matrix from a Matrix-Market file
  template <class T, class Prop, class Storage, class Allocator>
  void ReadMatrixMarket(string filename,
                        Matrix<T, Prop, Storage, Allocator>& A)
  {
    typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
    bool complex = ( sizeof(Tcplx)/sizeof(Treal) == 2);
    bool symmetric = IsSymmetricMatrix(A);
    
    ifstream input_stream(filename.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!input_stream.good())
      throw IOError("ReadMatrixMarket(string filename, Matrix& A)",
                    "Unable to read file \"" + filename + "\".");
#endif
    
    // first line
    string line, word1, word2, storage_type, value_type, matrix_type;
    getline(input_stream, line);
    
    istringstream header(line);
    header >> word1 >> word2 >> storage_type >> value_type >> matrix_type;
    
    if (storage_type != "coordinate")
      throw WrongArgument("ReadMatrixMarket(string filename, Matrix& A)",
                          "The storage should be coordinate in the file");
    
    if ( complex && (value_type != "complex") )
      throw WrongArgument("ReadMatrixMarket(string filename, Matrix& A)",
                          "The matrix should contain complex values");

    if ( !complex && (value_type != "real") )
      throw WrongArgument("ReadMatrixMarket(string filename, Matrix& A)",
                          "The matrix should contain real values");
    
    if ( symmetric && (matrix_type != "symmetric") )
      throw WrongArgument("ReadMatrixMarket(string filename, Matrix& A)",
                          "Problem of symmetry");

    if ( !symmetric && (matrix_type != "general") )
      throw WrongArgument("ReadMatrixMarket(string filename, Matrix& A)",
                          "Problem of symmetry");
    
    // skipping commentaries
    bool comment_line = true;
    while (comment_line)
      {
        getline(input_stream, line);
        if ( (line.size() > 0) && (line[0] != '%'))
          comment_line = false;
        
        if (input_stream.eof())
          comment_line = false;
      }
    
    // then reading m, n, nnz
    int m = 0, n = 0; long nnz = 0;
    istringstream size_stream(line);
    size_stream >> m >> n >> nnz;
    
    // then reading i, j, val
    Vector<int> row(nnz), col(nnz);
    Vector<Tcplx> val(nnz);
    ReadCoordinateMatrix(input_stream, row, col, val, complex);
    
    // closing stream
    input_stream.close();
    
    // for symmetric matrices, row and col are swapped
    // because Matrix Market is storing lower part whereas Seldon
    // stores upper part of the matrix
    A.Reallocate(m, n);
    if (symmetric)
      ConvertMatrix_from_Coordinates(col, row, val, A, 1);
    else
      ConvertMatrix_from_Coordinates(row, col, val, A, 1);
  }
  
  
  //! writing a matrix in Matrix-Market format
  template<class T, class Prop, class Storage, class Allocator>
  void WriteMatrixMarket(const Matrix<T, Prop, Storage, Allocator>& A,
                         const string& file_name)
  {
    typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type Tcplx;
    typedef typename ClassComplexType<T>::Treal Treal;
    bool complex = ( sizeof(Tcplx)/sizeof(Treal) == 2);
    bool symmetric = IsSymmetricMatrix(A);
    
    ofstream file_out(file_name.data());
    file_out.precision(15);
    file_out.setf(ios::scientific);
        
#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!file_out.good())
      throw IOError("WriteHarwellBoeing(const Matrix& A, string filename)",
                    "Unable to open file \"" + file_name + "\".");
#endif
    
    file_out << "%%MatrixMarket matrix coordinate ";
    if (complex)
      file_out << "complex ";
    else
      file_out << "real ";
    
    if (symmetric)
      file_out << "symmetric\n";
    else
      file_out << "general\n";
    
    // conversion to coordinate format (with 1-index)
    Vector<int> row, col;
    Vector<Tcplx> val;
    ConvertMatrix_to_Coordinates(A, row, col, val, 1, false);
    
    file_out << A.GetM() << "  " << A.GetN() << "  " << val.GetM() << '\n';
    
    if (symmetric)
      {
        // storing lower part for symmetric matrices
        int itmp;
        for (long i = 0; i < row.GetM(); i++)
          {
            if (row(i) < col(i))
              {
                itmp = row(i);
                row(i) = col(i);
                col(i) = itmp;
              }
          }
      }
    
    WriteCoordinateMatrix(file_out, row, col, val, complex);
    
    file_out.close();
  }
  
}


#define SELDON_FILE_MATRIX_SPARSE_IOMATRIXMARKET_CXX
#endif
