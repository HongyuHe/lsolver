// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_CXX

#include "Matrix_ArrayComplexSparse.hxx"

namespace Seldon
{


  //! Returns the number of non-zero elements (real part).
  /*!
    \return The number of non-zero elements for real part of matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  long Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetRealNonZeros() const
  {
    long nnz = 0;
    for (int i = 0; i < this->val_real_.GetM(); i++)
      nnz += this->val_real_(i).GetM();

    return nnz;
  }


  //! Returns the number of non-zero elements (imaginary part).
  /*!
    \return The number of non-zero elements for imaginary part of matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  long Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  GetImagNonZeros() const
  {
    int nnz = 0;
    for (int i = 0; i < this->val_imag_.GetM(); i++)
      nnz += this->val_imag_(i).GetM();

    return nnz;
  }


  //! returns size of matrix in bytes
  template<class T, class Prop, class Storage, class Allocator>
  size_t Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + sizeof(pointer)*this->val_real_.GetM();
    taille += sizeof(pointer)*this->val_imag_.GetM();
    for (int i = 0; i < val_real_.GetM(); i++)
      taille += val_real_(i).GetMemorySize();
    
    for (int i = 0; i < val_imag_.GetM(); i++)
      taille += val_imag_(i).GetMemorySize();
    
    return taille;
  }


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
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Reallocate(int i, int j)
  {
    // Clears previous entries.
    Clear();

    this->m_ = i;
    this->n_ = j;

    int n = Storage::GetFirst(i,j);
    val_real_.Reallocate(n);
    val_imag_.Reallocate(n);
  }


  //! Reallocates additional memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i number of rows.
    \param j number of columns.
    \note Data is kept.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  Resize(int i, int j)
  {
    int n = Storage::GetFirst(this->m_, this->n_);
    int new_n = Storage::GetFirst(i, j);
    if (n != new_n)
      {
	Vector<Vector<value_type, VectSparse, Allocator>, VectFull,
	  NewAlloc<Vector<value_type, VectSparse, Allocator> > > new_val_real;

	Vector<Vector<value_type, VectSparse, Allocator>, VectFull,
	  NewAlloc<Vector<value_type, VectSparse, Allocator> > > new_val_imag;

	new_val_real.Reallocate(new_n);
	new_val_imag.Reallocate(new_n);

	for (int k = 0 ; k < min(n, new_n) ; k++)
	  {
	    Swap(new_val_real(k), this->val_real_(k));
	    Swap(new_val_imag(k), this->val_imag_(k));
	  }

	val_real_.SetData(new_n, new_val_real.GetData());
	val_imag_.SetData(new_n, new_val_imag.GetData());
	new_val_real.Nullify();
	new_val_imag.Nullify();

      }

    this->m_ = i;
    this->n_ = j;
  }


  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets element (i, j) of matrix
  /*!
    \param[in] i row index.
    \param[in] j column index.
    \param[in] x A(i, j) = x
   */
  template <class T, class Prop, class Storage, class Allocator>
  void
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const entry_type& x)
  {
    if (real(x) != value_type(0))
      val_real_(Storage::GetFirst(i, j)).Get(Storage::GetSecond(i, j)) = real(x);
    else
      {
        if (val_real_(Storage::GetFirst(i, j))
	    (Storage::GetSecond(i, j)) != value_type(0))
          val_real_(Storage::GetFirst(i, j)).
	    Get(Storage::GetSecond(i, j)) = value_type(0);
      }
    
    if (imag(x) != value_type(0))
      val_imag_(Storage::GetFirst(i, j)).
	Get(Storage::GetSecond(i, j)) = imag(x);
    else
      {
        if (val_imag_(Storage::GetFirst(i, j))
	    (Storage::GetSecond(i, j)) != value_type(0))
          val_imag_(Storage::GetFirst(i, j)).
	    Get(Storage::GetSecond(i, j)) = value_type(0);
      }

  }
  
  
  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Print() const
  {
    if (Storage::GetFirst(1, 0) == 1)
      for (int i = 0; i < this->m_; i++)
	{
	  for (int j = 0; j < this->val_real_(i).GetM(); j++)
	    cout << (i+1) << " " << this->val_real_(i).Index(j)+1
		 << " " << this->val_real_(i).Value(j) << endl;

	  for (int j = 0; j < this->val_imag_(i).GetM(); j++)
	    cout << (i+1) << " " << this->val_imag_(i).Index(j)+1
		 << " (0, " << this->val_imag_(i).Value(j) << ")"<<endl;
	}
    else
      for (int i = 0; i < this->n_; i++)
	{
	  for (int j = 0; j < this->val_real_(i).GetM(); j++)
	    cout << this->val_real_(i).Index(j)+1 << " " << i+1
		 << " " << this->val_real_(i).Value(j) << endl;

	  for (int j = 0; j < this->val_imag_(i).GetM(); j++)
	    cout << this->val_imag_(i).Index(j)+1 << " " << i+1
		 << " (0, " << this->val_imag_(i).Value(j) << ")"<<endl;
	}
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
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArrayComplexSparse::Write(string FileName)",
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
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArrayComplexSparse::Write(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));

    for (int i = 0; i < val_real_.GetM(); i++)
      {
        val_real_(i).Write(FileStream);
        val_imag_(i).Write(FileStream);
      }
  }

  
  //! Writes the matrix in a file.
  /*! Stores the matrix in a file in ascii format. The entries are written in
    coordinate format (row column value). 1-index convention is used.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::WriteText(string FileName, bool cplx) const
  {
    ofstream FileStream; FileStream.precision(14);
    FileStream.open(FileName.c_str());
    
    // changing precision
    FileStream.precision(cout.precision());
    
#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArrayComplexSparse::WriteText(string FileName)",
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
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream, bool cplx) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_ArrayComplexSparse::"
                    "WriteText(ofstream& FileStream)",
		    "Stream is not ready.");
#endif

    // Conversion to coordinate format (1-index convention).
    const Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this);

    entry_type zero; int index = 1;
    WriteCoordinateMatrix(leaf_class, FileStream, zero, index, cplx);
  }

  
  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Read(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_ArrayComplexSparse::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    \param FileStream input stream
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream)
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

    val_real_.Reallocate(Storage::GetFirst(this->m_, this->n_));
    val_imag_.Reallocate(Storage::GetFirst(this->m_, this->n_));
    for (int i = 0; i < val_real_.GetM(); i++)
      {
        val_real_(i).Read(FileStream);
        val_imag_(i).Read(FileStream);
      }

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
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ReadText(string FileName, bool cplx)
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
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream, bool cplx)
  {
    Matrix<T, Prop, Storage, Allocator>& leaf_class =
      static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);
    
    entry_type zero; int index = 1;
    ReadCoordinateMatrix(leaf_class, FileStream, zero, index, -1, cplx);
  }
  
  
  //! Assembles the matrix.
  /*!
    All the row numbers are sorted.
    If same row numbers exist, values are added.
    \warning If you are using the methods AddInteraction/AddInteractions,
    you don't need to call that method.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Assemble()
  {
    for (int i = 0; i < Storage::GetFirst(this->m_, this->n_); i++)
      {
	val_real_(i).Assemble();
	val_imag_(i).Assemble();
      }
  }

  
  //! removes small entries
  template <class T, class Prop, class Storage, class Allocator> template<class T0>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>
  ::RemoveSmallEntry(const T0& epsilon)
  {
    for (int i = 0; i < val_real_.GetM(); i++)
      val_real_(i).RemoveSmallEntry(epsilon);

    for (int i = 0; i < val_imag_.GetM(); i++)
      val_imag_(i).RemoveSmallEntry(epsilon);
  }
  
  
  //! Matrix is initialized to the identity matrix.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::
  SetIdentity()
  {
    this->n_ = this->m_;
    for (int i = 0; i < this->m_; i++)
      {
        val_imag_(i).Clear();
	val_real_(i).Reallocate(1);
	val_real_(i).Index(0) = i;
	val_real_(i).Value(0) = value_type(1);
      }
  }


  //! Non-zero entries are set to 0 (but not removed).
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Zero()
  {
    for (int i = 0; i < Storage::GetFirst(this->m_, this->n_); i++)
      {
	val_real_(i).Zero();
	val_imag_(i).Zero();
      }
  }


  //! Non-zero entries are filled with values 0, 1, 2, 3 ...
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::Fill()
  {
    long value = 0;
    for (int i = 0; i < Storage::GetFirst(this->m_, this->n_); i++)
      {
	for (int j = 0; j < val_real_(i).GetM(); j++)
	  SetComplexReal(value++, val_real_(i).Value(j));

	for (int j = 0; j < val_imag_(i).GetM(); j++)
	  SetComplexReal(value++, val_imag_(i).Value(j));
      }
  }


  //! Non-zero entries are set to a given value x.
  /*!
    real non-zero entries are set to real(x)
    whereas imaginary non-zero entries are set to imag(x)
   */
  template <class T, class Prop, class Storage, class Allo> template<class T0>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allo>::
  Fill(const complex<T0>& x)
  {
    for (int i = 0; i < Storage::GetFirst(this->m_, this->n_); i++)
      {
	val_real_(i).Fill(real(x));
	val_imag_(i).Fill(imag(x));
      }
  }


  //! Non-zero entries are set to a given value x.
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>&
  Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::operator=
  (const complex<T0>& x)
  {
    this->Fill(x);
  }


  //! Non-zero entries take a random value.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_ArrayComplexSparse<T, Prop, Storage, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (int i = 0; i < Storage::GetFirst(this->m_, this->n_); i++)
      {
	val_real_(i).FillRand();
	val_imag_(i).FillRand();
      }
  }


  ////////////////////////////////////
  // MATRIX<ARRAY_COLCOMPLEXSPARSE> //
  ////////////////////////////////////


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<entry_type>& val, bool sorted)
  {
    for (int j = 0; j < nb; j++)
      AddInteraction(i, col(j), val(j));
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColComplexSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<entry_type>& val, bool sorted)
  {
    int nb_real = 0;
    int nb_imag = 0;
    IVect row_real(nb), row_imag(nb);
    Vector<value_type> val_real(nb), val_imag(nb);
    for (int j = 0; j < nb; j++)
      {
	if (real(val(j)) != value_type(0))
	  {
	    row_real(nb_real) = row(j);
	    val_real(nb_real) = real(val(j));
	    nb_real++;
	  }

	if (imag(val(j)) != value_type(0))
	  {
	    row_imag(nb_imag) = row(j);
	    val_imag(nb_imag) = imag(val(j));
	    nb_imag++;
	  }
      }

    this->val_real_(i).AddInteractionRow(nb_real, row_real, val_real, sorted);
    this->val_imag_(i).AddInteractionRow(nb_imag, row_imag, val_imag, sorted);
  }


  ////////////////////////////////////
  // MATRIX<ARRAY_ROWCOMPLEXSPARSE> //
  ////////////////////////////////////


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<entry_type>& val, bool sorted)
  {
    if (nb <= 0)
      return;

    int nb_real = 0;
    int nb_imag = 0;
    IVect col_real(nb), col_imag(nb);
    Vector<value_type> val_real(nb), val_imag(nb);
    for (int j = 0; j < nb; j++)
      {
	if (real(val(j)) != value_type(0))
	  {
	    col_real(nb_real) = col(j);
	    val_real(nb_real) = real(val(j));
	    nb_real++;
	  }

	if (imag(val(j)) != value_type(0))
	  {
	    col_imag(nb_imag) = col(j);
	    val_imag(nb_imag) = imag(val(j));
	    nb_imag++;
	  }
      }

    this->val_real_(i).AddInteractionRow(nb_real, col_real, val_real, sorted);
    this->val_imag_(i).AddInteractionRow(nb_imag, col_imag, val_imag, sorted);
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowComplexSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<entry_type>& val, bool sorted)
  {
    for (int j = 0; j < nb; j++)
      AddInteraction(row(j), i, val(j));
  }


  ///////////////////////////////////////
  // MATRIX<ARRAY_COLSYMCOMPLEXSPARSE> //
  ///////////////////////////////////////


  //! Sets element (i, j) of the matrix
  /*!
    \param i row index
    \param j column index
    \param x A(i, j) = x
   */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>
  ::Set(int i, int j, const entry_type& x)
  {
    if (i <= j)
      {
        if (real(x) != value_type(0))
          this->val_real_(j).Get(i) = real(x);
        else
          {
            if (this->val_real_(j)(i) != value_type(0))
              this->val_real_(j).Get(i) = value_type(0);
          }
        
        if (imag(x) != value_type(0))
          this->val_imag_(j).Get(i) = imag(x);
        else
          {
            if (this->val_imag_(j)(i) != value_type(0))
              this->val_imag_(j).Get(i) = value_type(0);
          }
      }
    else
      {
        if (real(x) != value_type(0))
          this->val_real_(i).Get(j) = real(x);
        else
          {
            if (this->val_real_(i)(j) != value_type(0))
              this->val_real_(i).Get(j) = value_type(0);
          }
        
        if (imag(x) != value_type(0))
          this->val_imag_(i).Get(j) = imag(x);
        else
          {
            if (this->val_imag_(i)(j) != value_type(0))
              this->val_imag_(i).Get(j) = value_type(0);
          }
      }
  }

      
  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  AddInteraction(int i, int j, const entry_type& val)
  {
    if (i <= j)
      {
	if (real(val) != value_type(0))
	  this->val_real_(j).AddInteraction(i, real(val));

	if (imag(val) != value_type(0))
	  this->val_imag_(j).AddInteraction(i, imag(val));
      }
  }


  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<entry_type>& val, bool sorted)
  {
    for (int j = 0; j < nb; j++)
      AddInteraction(i, col(j), val(j));
  }

  
  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<entry_type>& val, bool sorted)
  {
    int nb_real = 0;
    int nb_imag = 0;
    IVect row_real(nb), row_imag(nb);
    Vector<value_type> val_real(nb), val_imag(nb);
    for (int j = 0; j < nb; j++)
      if (row(j) <= i)
	{
	  if (real(val(j)) != value_type(0))
	    {
	      row_real(nb_real) = row(j);
	      val_real(nb_real) = real(val(j));
	      nb_real++;
	    }

	  if (imag(val(j)) != value_type(0))
	    {
	      row_imag(nb_imag) = row(j);
	      val_imag(nb_imag) = imag(val(j));
	      nb_imag++;
	    }
	}

    this->val_real_(i).AddInteractionRow(nb_real, row_real, val_real, sorted);
    this->val_imag_(i).AddInteractionRow(nb_imag, row_imag, val_imag, sorted);
  }


  ///////////////////////////////////////
  // MATRIX<ARRAY_ROWSYMCOMPLEXSPARSE> //
  ///////////////////////////////////////


  //! Sets element (i, j) of the matrix
  /*!
    \param i row index
    \param j column index
    \param x A(i, j) = x
   */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>
  ::Set(int i, int j, const entry_type& x)
  {
    if (i <= j)
      {
        if (real(x) != value_type(0))
          this->val_real_(i).Get(j) = real(x);
        else
          {
            if (this->val_real_(i)(j) != value_type(0))
              this->val_real_(i).Get(j) = value_type(0);
          }
        
        if (imag(x) != value_type(0))
          this->val_imag_(i).Get(j) = imag(x);
        else
          {
            if (this->val_imag_(i)(j) != value_type(0))
              this->val_imag_(i).Get(j) = value_type(0);
          }
      }
    else
      {
        if (real(x) != value_type(0))
          this->val_real_(j).Get(i) = real(x);
        else
          {
            if (this->val_real_(j)(i) != value_type(0))
              this->val_real_(j).Get(i) = value_type(0);
          }
        
        if (imag(x) != value_type(0))
          this->val_imag_(j).Get(i) = imag(x);
        else
          {
            if (this->val_imag_(j)(i) != value_type(0))
              this->val_imag_(j).Get(i) = value_type(0);
          }
      }
  }
  
  
  //! Adds a coefficient in the matrix.
  /*!
    \param[in] i row number.
    \param[in] j column number.
    \param[in] val coefficient to add.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  AddInteraction(int i, int j, const entry_type& val)
  {
    if (i <= j)
      {
	if (real(val) != value_type(0))
	  this->val_real_(i).AddInteraction(j, real(val));

	if (imag(val) != value_type(0))
	  this->val_imag_(i).AddInteraction(j, imag(val));
      }
  }
  

  //! Adds coefficients in a row.
  /*!
    \param[in] i row number.
    \param[in] nb number of coefficients to add.
    \param[in] col column numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  AddInteractionRow(int i, int nb, const IVect& col,
		    const Vector<entry_type>& val, bool sorted)
  {
    int nb_real = 0;
    int nb_imag = 0;
    IVect col_real(nb), col_imag(nb);
    Vector<value_type> val_real(nb), val_imag(nb);
    for (int j = 0; j < nb; j++)
      if (i <= col(j))
	{
	  if (real(val(j)) != value_type(0))
	    {
	      col_real(nb_real) = col(j);
	      val_real(nb_real) = real(val(j));
	      nb_real++;
	    }

	  if (imag(val(j)) != value_type(0))
	    {
	      col_imag(nb_imag) = col(j);
	      val_imag(nb_imag) = imag(val(j));
	      nb_imag++;
	    }
	}

    this->val_real_(i).AddInteractionRow(nb_real, col_real, val_real, sorted);
    this->val_imag_(i).AddInteractionRow(nb_imag, col_imag, val_imag, sorted);
  }


  //! Adds coefficients in a column.
  /*!
    \param[in] i column number.
    \param[in] nb number of coefficients to add.
    \param[in] row row numbers of coefficients.
    \param[in] val values of coefficients.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>::
  AddInteractionColumn(int i, int nb, const IVect& row,
		       const Vector<entry_type>& val, bool sorted)
  {
    for (int j = 0; j < nb; j++)
      AddInteraction(row(j), i, val(j));
  }

} // namespace Seldon

#define SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_CXX
#endif
