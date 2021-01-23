// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_ARRAY_SPARSE_CXX

#include "Matrix_ArraySparse.hxx"

namespace Seldon
{


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    \warning Data is lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Reallocate(int i, int j)
  {
    // Clears previous entries.
    Clear();

    this->m_ = i;
    this->n_ = j;

    int n = Storage::GetFirst(i, j);
    val_.Reallocate(n);
  }

  
  //! Reallocates additional memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    Data is kept
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Resize(int i, int j)
  {
    int n = Storage::GetFirst(this->m_, this->n_);
    int new_n = Storage::GetFirst(i, j);
    if (n != new_n)
      {
	Vector<Vector<T, VectSparse, Allocator>, VectFull,
	  NewAlloc<Vector<T, VectSparse, Allocator> > > new_val;

	new_val.Reallocate(new_n);

	for (int k = 0 ; k < min(n, new_n) ; k++)
	  Swap(new_val(k), this->val_(k));

	val_.SetData(new_n, new_val.GetData());
	new_val.Nullify();

      }

    this->m_ = i;
    this->n_ = j;
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of non-zero entries.
  /*!
    \return The number of non-zero entries.
  */
  template <class T, class Prop, class Storage, class Allocator>
  long Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetNonZeros()
    const
  {
    long nnz = 0;
    for (int i = 0; i < this->val_.GetM(); i++)
      nnz += this->val_(i).GetM();

    return nnz;
  }


  //! returns size of matrix in bytes
  template<class T, class Prop, class Storage, class Allocator>
  size_t Matrix_ArraySparse<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + sizeof(pointer)*this->val_.GetM();
    for (int i = 0; i < this->val_.GetM(); i++)
      taille += this->val_(i).GetMemorySize();
    
    return taille;
  }
  

  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Print() const
  {
    if (Storage::GetFirst(1, 0) == 1)
      for (int i = 0; i < this->m_; i++)
	{
	  for (int j = 0; j < this->val_(i).GetM(); j++)
	    cout << (i+1) << " " << this->val_(i).Index(j)+1
		 << " " << this->val_(i).Value(j) << endl;
	}
    else
      for (int i = 0; i < this->n_; i++)
	{
	  for (int j = 0; j < this->val_(i).GetM(); j++)
	    cout << this->val_(i).Index(j)+1 << " " << i+1
		 << " " << this->val_(i).Value(j) << endl;
	}
  }


  //! Assembles the matrix
  /*!
    All the column/row numbers are sorted.
    If same column/row numbers exist, values are added.
    \warning If you are using the methods AddInteraction,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Assemble()
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Assemble();
  }


  //! Removes small coefficients from entries.
  /*!
    \param[in] epsilon entries whose values are below epsilon are removed.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  RemoveSmallEntry(const T0& epsilon)
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).RemoveSmallEntry(epsilon);
  }


  //! Matrix is initialized to the identity matrix.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::SetIdentity()
  {
    T one; SetComplexOne(one);
    this->n_ = this->m_;
    for (int i = 0; i < this->m_; i++)
      {
	val_(i).Reallocate(1);
	val_(i).Index(0) = i;
	val_(i).Value(0) = one;
      }
  }


  //! Non-zero entries are set to 0 (but not removed).
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Zero()
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Zero();
  }


  //! Non-zero entries are filled with values 0, 1, 2, 3 ...
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::Fill()
  {
    int value = 0;
    for (int i = 0; i < val_.GetM(); i++)
      for (int j = 0; j < val_(i).GetM(); j++)
	{
          SetComplexReal(value, val_(i).Value(j));
          value++;
        }
  }


  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allo> template<class T0>
  void Matrix_ArraySparse<T, Prop, Storage, Allo>::Fill(const T0& x)
  {
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Fill(x);
  }


  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_ArraySparse<T, Prop, Storage, Allocator>&
  Matrix_ArraySparse<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);
  }


  //! Non-zero entries take a random value.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).FillRand();
  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Write(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));

    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Write(FileStream);
  }


  //! Writes the matrix in a file.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row column value). 1-index convention is used.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  WriteText(string FileName, bool cplx) const
  {
    ofstream FileStream; FileStream.precision(14);
    FileStream.open(FileName.c_str());

    // changing precision
    FileStream.precision(cout.precision());
    
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream, cplx);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row column value). 1-index convention is used.
    \param FileStream output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  WriteText(ostream& FileStream, bool cplx) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    // conversion in coordinate format (1-index convention)
    const Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this);

    T zero; int index = 1;
    WriteCoordinateMatrix(leaf_class, FileStream, zero, index, cplx);
  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  Read(istream& FileStream)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Read(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.read(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		    sizeof(int));
    FileStream.read(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		    sizeof(int));

    val_.Reallocate(Storage::GetFirst(this->m_, this->n_));
    for (int i = 0; i < val_.GetM(); i++)
      val_(i).Read(FileStream);

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_ArraySparse::Read(istream& FileStream)",
                    string("Input operation failed.")
		    + string(" The input file may have been removed")
		    + " or may not contain enough data.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads the matrix from a file in text format.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ReadText(string FileName, bool cplx)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArraySparse::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream, cplx);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix from a stream in text format.
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArraySparse<T, Prop, Storage, Allocator>::
  ReadText(istream& FileStream, bool cplx)
  {
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);
    
    T zero; int index = 1;
    ReadCoordinateMatrix(leaf_class, FileStream, zero, index, -1, cplx);
  }


  /////////////////////////////
  // MATRIX<ARRAY_COLSPARSE> //
  /////////////////////////////


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_,
		    bool already_sorted)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val, already_sorted);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_,
                       bool already_sorted)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val, already_sorted);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionRow(int i, int nb, const Vector<int>& col,
		    const Vector<T>& val, bool already_sorted)
  {
    for (int j = 0; j < nb; j++)
      this->val_(col(j)).AddInteraction(i, val(j));
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const Vector<int>& row,
		       const Vector<T>& val, bool already_sorted)
  {
    this->val_(i).AddInteractionRow(nb, row, val, already_sorted);
  }


  /////////////////////////////
  // MATRIX<ARRAY_ROWSPARSE> //
  /////////////////////////////


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_,
                    bool already_sorted)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val, already_sorted);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_,
		       bool already_sorted)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val, already_sorted);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionRow(int i, int nb, const Vector<int>& col,
		    const Vector<T>& val, bool already_sorted)
  {
    this->val_(i).AddInteractionRow(nb, col, val, already_sorted);
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const Vector<int>& row,
		       const Vector<T>& val, bool already_sorted)
  {
    for (int j = 0; j < nb; j++)
      this->val_(row(j)).AddInteraction(i, val(j));
  }


  ////////////////////////////////
  // MATRIX<ARRAY_COLSYMSPARSE> //
  ////////////////////////////////


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_,
		    bool already_sorted)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val, already_sorted);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_,
                       bool already_sorted)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val, already_sorted);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, const Vector<int>& col,
		    const Vector<T>& val, bool already_sorted)
  {
    for (int j = 0; j < nb; j++)
      if (i <= col(j))
	this->val_(col(j)).AddInteraction(i, val(j));
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const Vector<int>& row,
		       const Vector<T>& val, bool already_sorted)
  {
    IVect new_row(nb);
    Vector<T> new_val(nb);
    nb = 0;
    for (int j = 0; j < new_row.GetM(); j++)
      if (row(j) <= i)
	{
	  new_row(nb) = row(j);
	  new_val(nb) = val(j); nb++;
	}

    this->val_(i).AddInteractionRow(nb, new_row, new_val, already_sorted);
  }

  
  ////////////////////////////////
  // MATRIX<ARRAY_ROWSYMSPARSE> //
  ////////////////////////////////


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col_ column numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, int* col_, T* value_,
                    bool already_sorted)
  {
    IVect col;
    col.SetData(nb, col_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionRow(i, nb, col, val, already_sorted);
    col.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row_ row numbers of coefficients.
    \param[in] value_ values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, int* row_, T* value_,
		       bool already_sorted)
  {
    IVect row;
    row.SetData(nb, row_);
    Vector<T> val;
    val.SetData(nb, value_);
    AddInteractionColumn(i, nb, row, val, already_sorted);
    row.Nullify();
    val.Nullify();
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionRow(int i, int nb, const Vector<int>& col,
		    const Vector<T>& val, bool already_sorted)
  {
    IVect new_col(nb);
    Vector<T> new_val(nb);
    nb = 0;
    for (int j = 0; j < new_col.GetM(); j++)
      if (i <= col(j))
	{
	  new_col(nb) = col(j);
	  new_val(nb) = val(j); nb++;
	}

    this->val_(i).AddInteractionRow(nb, new_col, new_val, already_sorted);
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const Vector<int>& row,
		       const Vector<T>& val, bool already_sorted)
  {
    for (int j = 0; j < nb; j++)
      if (row(j) <= i)
        this->val_(row(j)).AddInteraction(i, val(j));
  }
  
} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_SPARSE_CXX
#endif
