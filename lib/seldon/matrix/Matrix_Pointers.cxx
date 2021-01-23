// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_POINTERS_CXX

#include "Matrix_Pointers.hxx"

namespace Seldon
{
  

  //! Main constructor.
  /*! Builds a i x j full matrix.
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Matrix_Pointers(int i, int j): Matrix_Base<T, Allocator>(i, j)
  {

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>( calloc(Storage::GetFirst(i, j),
						 sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
      }
    if (me_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Pointers::Matrix_Pointers(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->data_ = Allocator::allocate(long(i) * long(j), this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(me_);
	me_ = NULL;
	this->data_ = NULL;
      }
    if (this->data_ == NULL && i != 0 && j != 0)
      throw NoMemory("Matrix_Pointers::Matrix_Pointers(int, int)",
		     string("Unable to allocate memory for a matrix of size ")
		     + to_str(static_cast<long int>(i)
			      * static_cast<long int>(j)
			      * static_cast<long int>(sizeof(T)))
		     + " bytes (" + to_str(i) + " x " + to_str(j)
		     + " elements).");
#endif

    pointer ptr = this->data_;
    long lgth = Storage::GetSecond(i, j);
    for (int k = 0; k < Storage::GetFirst(i, j); k++, ptr += lgth)
      me_[k] = ptr;

  }


  //! Copy constructor.
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Matrix_Pointers(const Matrix_Pointers<T, Prop, Storage, Allocator>& A):
    Matrix_Base<T, Allocator>(A)
  {
    this->m_ = 0;
    this->n_ = 0;
    this->data_ = NULL;
    this->me_ = NULL;

    this->Copy(A);
  }


  //! Clears the matrix.
  /*!
    Destructs the matrix.
    \warning On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Clear()
  {
#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (this->data_ != NULL)
	  {
	    Allocator::deallocate(this->data_, long(this->m_) * long(this->n_));
	    this->data_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (me_ != NULL)
	  {
	    free(me_);
	    me_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
      }
#endif
    
    this->m_ = 0;
    this->n_ = 0;
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Reallocates memory to resize the matrix.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning Depending on your allocator, data may be lost.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Reallocate(int i, int j)
  {

    if (i != this->m_ || j != this->n_)
      {
	this->m_ = i;
	this->n_ = j;

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif

	    me_ = reinterpret_cast<pointer*>( realloc(me_,
						      Storage::GetFirst(i, j)
						      * sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (me_ == NULL && i != 0 && j != 0)
	  throw NoMemory("Matrix_Pointers::Reallocate(int, int)",
			 string("Unable to reallocate memory for")
			 + " a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " elements).");
#endif

#ifdef SELDON_CHECK_MEMORY
	try
	  {
#endif
            long taille = long(i)*long(j);
	    this->data_ =
	      reinterpret_cast<pointer>(Allocator::
					reallocate(this->data_, taille,
						   this));

#ifdef SELDON_CHECK_MEMORY
	  }
	catch (...)
	  {
	    this->m_ = 0;
	    this->n_ = 0;
	    free(me_);
	    me_ = NULL;
	    this->data_ = NULL;
	  }
	if (this->data_ == NULL && i != 0 && j != 0)
	  throw NoMemory("Matrix_Pointers::Reallocate(int, int)",
			 string("Unable to reallocate memory")
			 + " for a matrix of size "
			 + to_str(static_cast<long int>(i)
				  * static_cast<long int>(j)
				  * static_cast<long int>(sizeof(T)))
			 + " bytes (" + to_str(i) + " x " + to_str(j)
			 + " elements).");
#endif

	pointer ptr = this->data_;
	long lgth = Storage::GetSecond(i, j);
	for (int k = 0; k < Storage::GetFirst(i, j); k++, ptr += lgth)
	  me_[k] = ptr;
      }
  }


  //! Changes the size of the matrix and sets its data array
  //! (low level method).
  /*!
    The matrix is first cleared (memory is freed). The matrix is then resized
    to a i x j matrix, and the data array of the matrix is set to 'data'.
    'data' elements are not duplicated: the new data array of the matrix is
    the 'data' array. It is useful to create a matrix from pre-existing data.
    \param i new number of rows.
    \param j new number of columns.
    \param data new array storing elements.
    \warning 'data' has to be used carefully outside the object.
    Unless you use 'Nullify', 'data' will be freed by the destructor,
    which means that 'data' must have been allocated carefully. The matrix
    allocator should be compatible.
    \note This method should only be used by advanced users.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::SetData(int i, int j,
	    typename Matrix_Pointers<T, Prop, Storage, Allocator>
	    ::pointer data)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = j;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	me_ = reinterpret_cast<pointer*>( calloc(Storage::GetFirst(i, j),
						 sizeof(pointer)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
	this->data_ = NULL;
	return;
      }
    if (me_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	this->data_ = NULL;
	return;
      }
#endif

    this->data_ = data;

    pointer ptr = this->data_;
    long lgth = Storage::GetSecond(i, j);
    for (int k = 0; k < Storage::GetFirst(i, j); k++, ptr += lgth)
      me_[k] = ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
    \warning Memory is not released except for me_.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Nullify()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	if (me_ != NULL)
	  {
	    free(me_);
	    me_ = NULL;
	  }

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	me_ = NULL;
      }
#endif

    this->data_ = NULL;
  }

  
  //! Reallocates memory to resize the matrix and keeps previous entries.
  /*!
    On exit, the matrix is a i x j matrix.
    \param i new number of rows.
    \param j new number of columns.
    \warning The previous entries are kept, extra-entries may not be
    initialized (depending of the allocator).
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Resize(int i, int j)
  {

    if (i == this->m_ && j == this->n_)
      return;

    // Storing the old values of the matrix.
    int iold = Storage::GetFirst(this->m_, this->n_);
    int jold = Storage::GetSecond(this->m_, this->n_);
    Vector<value_type, VectFull, Allocator> xold(this->GetDataSize());
    for (long k = 0; k < this->GetDataSize(); k++)
      xold(k) = this->data_[k];

    // Reallocation.
    int inew = Storage::GetFirst(i, j);
    int jnew = Storage::GetSecond(i, j);
    this->Reallocate(i,j);

    // Filling the matrix with its old values.
    int imin = min(iold, inew), jmin = min(jold, jnew);
    for (int k = 0; k < imin; k++)
      for (int l = 0; l < jmin; l++)
	this->data_[k*jnew+l] = xold(l+jold*k);
  }

  
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Sets all elements to zero.
  /*!
    \warning It fills the memory with zeros. If the matrix stores complex
    structures, use 'Fill' instead.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Zero()
  {
    Allocator::memoryset(this->data_, char(0),
			 this->GetDataSize() * sizeof(value_type));
  }


  //! Sets the current matrix to the identity.
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::SetIdentity()
  {
    T zero, one;
    SetComplexZero(zero);
    SetComplexOne(one);
    
    Fill(zero);
    
    for (int i = 0; i < min(this->m_, this->n_); i++)
      (*this)(i,i) = one;
  }


  //! Fills the matrix the matrix with 0, 1, 2, ...
  /*!
    On exit, the matrix is filled with 0, 1, 2, 3, ... The order of
    those numbers depends on the storage.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Fill()
  {
    for (long i = 0; i < this->GetDataSize(); i++)
      SetComplexReal(i,  this->data_[i]);
  }


  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator> 
  template <class T0>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    T x_; SetComplexReal(x, x_);
    for (long i = 0; i < this->GetDataSize(); i++)
      this->data_[i] = x_;
  }

  //! Fills the matrix with a given value.
  /*!
    \param x the value to fill the matrix with.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class T0>
  Matrix_Pointers<T, Prop, Storage, Allocator>&
  Matrix_Pointers<T, Prop, Storage, Allocator>::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Fills the matrix randomly.
  /*!
    \note The random generator is very basic.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::FillRand()
  {
#ifndef SELDON_WITHOUT_REINIT_RANDOM
    srand(time(NULL));
#endif
    for (long i = 0; i < this->GetDataSize(); i++)
      SetComplexReal(rand(), this->data_[i]);
  }


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


  //! Displays a sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its upper-left corner (a, b)
    and its bottom-right corner (m, n). So, elements with indices
    in [a, m] x [b, n] are displayed on the standard output,
    in text format. Each row is displayed on a single line and
    elements of a row are delimited by tabulations.
    \param a row index of the upper-left corner.
    \param b column index of the upper-left corner.
    \param m row index of the bottom-right corner.
    \param n column index of the bottom-right corner.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Print(int a, int b,
							   int m, int n) const
  {
    for (int i = a; i < min(this->m_, a+m); i++)
      {
	for (int j = b; j < min(this->n_, b+n); j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


  //! Displays a square sub-matrix on the standard output.
  /*!
    The sub-matrix is defined by its bottom-right corner (l, l).
    So, elements with indices in [0, 0] x [l, l] are displayed
    on the standard output, in text format. Each row is displayed
    on a single line and elements of a row are delimited
    by tabulations.
    \param l dimension of the square matrix to be displayed.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::Print(int l) const
  {
    Print(0, 0, l, l);
  }


  /**************************
   * INPUT/OUTPUT FUNCTIONS *
   **************************/

#ifdef SELDON_WITH_HDF5
  //! Writes the matrix in a HDF5 file.
  /*!
    All elements of the matrix are stored in HDF5 format.
    \param FileName file name.
    \param group_name name of the group the matrix must be stored in.
    \param dataset_name name of the dataset the matrix must be stored in.
    \warning This method is not defined yet!
  */
 template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
 ::WriteHDF5(string FileName, string group_name, string dataset_name) const
  {
    throw IOError("Matrix_Pointers::WriteHDF5(string FileName)",
                  string("Unable to write matrix in \"") + FileName + "\".");
  }
#endif


  //! Appends the matrix in a file.
  /*!  
    Appends the matrix in a file in binary format. The matrix elements are
    appended in the same order as in memory (e.g. row-major storage).
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Append(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary | ios::app);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, false);

    FileStream.close();
  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileName output file name.
    \param with_size if set to 'false', the dimensions of the matrix are not
    saved.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Write(string FileName, bool with_size) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str(), ofstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream, with_size);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in binary format.
    The number of rows (integer) and the number of columns (integer)
    are written, and matrix elements are then written in the same order
    as in memory (e.g. row-major storage).
    \param FileStream output stream.
    \param with_size if set to 'false', the dimensions of the matrix are not
    saved.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream, bool with_size) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Write(ostream& FileStream)",
		    "The stream is not ready.");
#endif

    if (with_size)
      {
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
			 sizeof(int));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
			 sizeof(int));
      }

    FileStream.write(reinterpret_cast<char*>(this->data_),
		     long(this->m_) * long(this->n_) * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Write(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Writes the matrix in a file.
  /*!
    Stores the matrix in a file in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileName output file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }


  //! Writes the matrix to an output stream.
  /*!
    Writes the matrix to an output stream in text format.
    Only matrix elements are written (not dimensions).
    Each row is written on a single line and elements of
    a row are delimited by tabulations.
    \param FileStream output stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {

#ifdef SELDON_CHECK_IO
    // Checks if the file is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::WriteText(ostream& FileStream)",
                    "The stream is not ready.");
#endif

    int i, j;
    for (i = 0; i < this->GetM(); i++)
      {
	for (j = 0; j < this->GetN(); j++)
	  FileStream << (*this)(i, j) << '\t';
	FileStream << endl;
      }

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::WriteText(ostream& FileStream)",
                    "Output operation failed.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix stored in binary format in a file.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileName input file name.
    \param with_size if set to 'false', the dimensions of the matrix are not
    available in the file. In this case, the current dimensions (M, N) of the
    matrix are unchanged, and MxN elements are read in the file.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Read(string FileName, bool with_size)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str(), ifstream::binary);

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::Read(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Read(FileStream, with_size);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in binary format from an input stream.
    The number of rows (integer) and the number of columns (integer)
    are read, and matrix elements are then read in the same order
    as it should be in memory (e.g. row-major storage).
    \param FileStream input stream.
    \param with_size if set to 'false', the dimensions of the matrix are not
    available in the stream. In this case, the current dimensions (M, N) of
    the matrix are unchanged, and MxN elements are read in the stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::Read(istream& FileStream, bool with_size)
  {

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(istream& FileStream)",
                    "The stream is not ready.");
#endif

    if (with_size)
      {
        int new_m, new_n;
        FileStream.read(reinterpret_cast<char*>(&new_m), sizeof(int));
        FileStream.read(reinterpret_cast<char*>(&new_n), sizeof(int));
        this->Reallocate(new_m, new_n);
      }

    FileStream.read(reinterpret_cast<char*>(this->data_),
		    long(this->GetM()) * long(this->GetN())
                    * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was read.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::Read(istream& FileStream)",
                    "Input operation failed.");
#endif

  }


  //! Reads the matrix from a file.
  /*!
    Reads a matrix from a file in text format.
    \param FileName input file name.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>::ReadText(string FileName)
  {
    ifstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Pointers::ReadText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->ReadText(FileStream);

    FileStream.close();
  }


  //! Reads the matrix from an input stream.
  /*!
    Reads a matrix in text format from an input stream.
    \param FileStream input stream.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Pointers<T, Prop, Storage, Allocator>
  ::ReadText(istream& FileStream)
  {
    // Clears the matrix.
    Clear();

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix_Pointers::ReadText(istream& FileStream)",
                    "The stream is not ready.");
#endif

    // Reads the first line.
    string line;
    getline(FileStream, line);
    if (FileStream.fail())
      // Is the file empty?
      return;

    // Converts the first line into a vector.
    istringstream line_stream(line);
    Vector<T> first_row;
    first_row.ReadText(line_stream);

    // Now reads all other rows, and puts them in a single vector.
    Vector<T> other_row;
    other_row.ReadText(FileStream);
    
    // Number of rows and columns.
    int n = first_row.GetM();
    int m = 1 + other_row.GetM() / n;

#ifdef SELDON_CHECK_IO
    // Checks that enough elements were read.
    if (other_row.GetM() != (m - 1) * n)
      throw IOError("Matrix_Pointers::ReadText(istream& FileStream)",
                    "Not all rows have the same number of columns.");
#endif

    this->Reallocate(m, n);
    // Fills the matrix.
    for (int j = 0; j < n; j++)
      this->Val(0, j) = first_row(j);

    int k = 0;
    for (int i = 1; i < m; i++)
      for (int j = 0; j < n; j++)
	this->Val(i, j) = other_row(k++);
  }


  //////////////////////
  // MATRIX<COLMAJOR> //
  //////////////////////


  /*****************
   * OTHER METHODS *
   *****************/


  //! Writes one column of matrix in a file.
  /*!
    Stores one column of matrix in a file in binary format.
    The column indexed by \a col is written.
    \param FileName output file name.
    \param col column index.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ColMajor, Allocator>
  ::WriteColumn(string FileName, int col) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix::WriteColumn(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteColumn(FileStream, col);

    FileStream.close();
  }


  //! Writes one column of matrix in an output stream.
  /*!
    Stores one column of matrix in an output stream in binary format.
    The column indexed by \a col is written.
    \param FileStream output stream.
    \param col column index.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, ColMajor, Allocator>
  ::WriteColumn(ostream& FileStream, int col) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (col < 0 || col >= this->n_)
      throw WrongCol("Matrix::WriteColumn(ostream& FileStream, int col)",
		     string("Index should be in [0, ")
		     + to_str(this->n_-1) + "], but is equal to "
		     + to_str(col) + ".");
#endif

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix::WriteColumn(ostream& FileStream, int col)",
		    "The stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(this->me_[col]),
		     this->m_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix::WriteColumn(ostream& FileStream, int col)",
                    "Output operation failed.");
#endif

  }


  //////////////////////
  // MATRIX<ROWMAJOR> //
  //////////////////////


  /*****************
   * OTHER METHODS *
   *****************/


  //! Writes one row of matrix in a file.
  /*!
    Stores one row of matrix in a file in binary format.
    The row indexed by \a row is written.
    \param FileName output file name.
    \param row row index.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, RowMajor, Allocator>
  ::WriteRow(string FileName, int row) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix::WriteRow(string FileName, int row)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteRow(FileStream, row);

    FileStream.close();
  }


  //! Writes one row of matrix in an output stream.
  /*!
    Stores one row of matrix in an output stream in binary format.
    The row indexed by \a row is written.
    \param FileStream output stream.
    \param row row index.
  */
  template <class T, class Prop, class Allocator>
  void Matrix<T, Prop, RowMajor, Allocator>
  ::WriteRow(ostream& FileStream, int row) const
  {
#ifdef SELDON_CHECK_BOUNDS
    if (row < 0 || row >= this->m_)
      throw WrongRow("Matrix::WriteRow(ostream& FileStream, int row)",
		     string("Index should be in [0, ")
		     + to_str(this->m_-1) + "], but is equal to "
		     + to_str(row) + ".");
#endif

#ifdef SELDON_CHECK_IO
    // Checks if the stream is ready.
    if (!FileStream.good())
      throw IOError("Matrix::WriteRow(ostream& FileStream, int rowm)",
		    "The stream is not ready.");
#endif

    FileStream.write(reinterpret_cast<char*>(this->me_[row]),
		     this->n_ * sizeof(value_type));

#ifdef SELDON_CHECK_IO
    // Checks if data was written.
    if (!FileStream.good())
      throw IOError("Matrix::WriteRow(ostream& FileStream, int row)",
                    "Output operation failed.");
#endif

  }
  
  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_POINTERS_CXX
#endif
