// Copyright (C) 2014 INRIA
// Author(s): Marc Duruflé
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

#ifndef SELDON_FILE_BAND_MATRIX_CXX

#include "BandMatrixInline.cxx"

namespace Seldon
{

  /***************
   * Matrix_Band *
   ***************/
  
  
  //! default constructor
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_Band<T, Prop, Storage, Allocator>::Matrix_Band()
  {
    kl_ = 0;
    ku_ = 0;
    this->m_ = 0;
    this->n_ = 0;    
  }
  
  
  //! clears the matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::Clear()
  {
    kl_ = 0;
    ku_ = 0;
    this->m_ = 0;
    this->n_ = 0;
    data_.Clear();
  }
  
  
  //! changes the size of the matrix, previous entries are lost
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::Reallocate(int m, int n, int kl, int ku)
  {
    this->m_ = m;
    this->n_ = n;
    kl_ = kl;
    ku_ = ku;
    data_.Reallocate(2*kl_+ku_+1, this->n_);
    data_.Fill(0);
  }
  

  //! sets A(i, j) = A(i, j) + val
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    int k = kl_ + ku_ + i - j;
    if ((k >= kl_) && (k <= 2*kl_+ku_))
      data_(k, j) += val;
    else
      {
        cout << "Matrix not compatible " << endl;
	abort();
      }
  }
  
  
  //! adds severals values on a single row of the matrix  
  /*!
    \param[in] i row on which values are added
    \param[in] n number of values to add
    \param[in] num column numbers
    \param[in] val values
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::
  AddInteractionRow(int i, int n, const IVect& num, const Vector<T>& val,
		    bool sorted)
  {
    for (int j = 0; j < n; j++)
      AddInteraction(i, num(j), val(j));
  }
  
  
  //! clears row i, fills it with 0
  /*!
    \param[in] i row to clear
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::ClearRow(int i)
  {
    T zero; SetComplexZero(zero);
    for (int k = max(-i, -this->kl_); k <= min(this->ku_, this->n_-1-i); k++)
      {
        int j = i+k;
        this->Get(i, j) = zero;
      }
  }
  
  
  //! returns entry (i, j) of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  const T Matrix_Band<T, Prop, Storage, Allocator>
  ::operator()(int i, int j) const
  {
    T zero; SetComplexZero(zero);
    int k = kl_ + ku_ + i - j;
    if ((k >= kl_) && (k <= 2*kl_+ku_))
      return data_(k, j);
    
    return zero;
  }  
  
  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  T& Matrix_Band<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    int k = kl_ + ku_ + i - j;
    if ((k >= kl_) && (k <= 2*kl_+ku_))
      return data_(k, j);
    else
      {
	cout << "Element not accessible " << endl;
	abort();
      }
  }
  
  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  const T& Matrix_Band<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    int k = kl_ + ku_ + i - j;
    if ((k >= kl_) && (k <= 2*kl_+ku_))
      return data_(k, j);
    else
      {
	cout << "Element not accessible " << endl;
	abort();
      }
  }


  //! sets the matrix to the identity matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::SetIdentity()
  {
    T zero, one; SetComplexZero(zero);
    SetComplexOne(one);
    
    data_.Fill(zero);
    for (int i = 0; i < this->n_; i++)
      data_(kl_+ku_, i) = one;
  }
  
  
  //! sets all non-zero entries to a same value
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void Matrix_Band<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    data_.Fill(x);

    T zero; SetComplexZero(zero);
    for (int i = 0; i < this->n_; i++)
      for (int j = 0; j < kl_; j++)
	data_(j, i) = zero;
  }
  
  
  //! sets all non-zero entries to random values
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::FillRand()
  {
    data_.FillRand();

    T zero; SetComplexZero(zero);
    for (int i = 0; i < this->n_; i++)
      for (int j = 0; j < kl_; j++)
	data_(j, i) = zero;
  }
  
  
  //! conversion from a CSR matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::Copy(const Matrix<T, General, ArrayRowSparse>& A)
  {
    Clear();

    // finding kl and ku
    int kl = 0, ku = 0;
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        {
          int col = A.Index(i, j);
          if (col < i)
            kl = max(kl, i-col);
          else if (col > i)
            ku = max(ku, col-i);
        }
    
    Reallocate(A.GetM(), A.GetN(), kl, ku);
    T zero; SetComplexZero(zero);
    Fill(zero);
    
    // filling the matrix with values of A
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        AddInteraction(i, A.Index(i, j), A.Value(i, j));
    
  }
    
  
  //! performs LU factorisation without pivoting
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::Factorize()
  {
    // position of the main diagonal
    int d = kl_ + ku_;
    T pivot, one;
    SetComplexOne(one);
    for (int j = 0; j < this->n_; j++)
      {
        // replacing diagonal element by its inverse
        data_(d, j) = one/data_(d, j);
        for (int p = 1; p <= min(this->m_ - 1 - j, kl_); p++)
          {
            pivot = data_(d+p, j)*data_(d, j);
            data_(d+p, j) = pivot;
            for (int k = 1; k <= min(this->n_ - 1 - j, ku_); k++)
              data_(d+p-k, j+k) -= pivot*data_(d-k, j+k);
          }
      }
  }
  

  //! performs LU factorisation with pivoting
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::Factorize(IVect& ipivot)
  {
    ipivot.Reallocate(this->m_);
    // position of the main diagonal
    int d = kl_ + ku_;
    T pivot, one, tmp;
    typename ClassComplexType<T>::Treal amax;
    SetComplexOne(one);
    for (int j = 0; j < this->n_; j++)
      {
	// searching pivot
	int pmin = min(this->m_ - 1 - j, kl_);
	amax = abs(data_(d, j));
	int p0 = 0;
	for (int p = 1; p <= pmin; p++)
	  if (abs(data_(d+p, j)) > amax)
	    {
	      p0 = p;
	      amax = abs(data_(d+p, j));
	    }
	
	ipivot(j) = j+p0;
	if (p0 > 0)
	  {
	    // interchanging row j and j + p0
	    for (int j2 = j; j2 <= min(this->m_-1, kl_ + ku_ + j); j2++)
	      {
		int k = kl_ + ku_ + j - j2;
		tmp = data_(k, j2);
		data_(k, j2) = data_(k+p0, j2);
		data_(k+p0, j2) = tmp;
	      }	    
	  }
	
	// replacing diagonal element by its inverse
        data_(d, j) = one/data_(d, j);
	
	// performing Gauss elimination on column j
	for (int p = 1; p <= pmin; p++)
	  {
	    pivot = data_(d+p, j)*data_(d, j);
            data_(d+p, j) = pivot;
	    
            for (int k = 1; k <= min(this->n_ - 1 - j, kl_ + ku_); k++)
	      data_(d+p-k, j+k) -= pivot*data_(d-k, j+k);
	  }
      }
  }
  
  
  //! performs the operation *this = *this + alpha A
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0> void Matrix_Band<T, Prop, Storage, Allocator>
  ::Add_(const T0& alpha, const Matrix<T, General, BandedCol, Allocator>& A)
  {
    int kl = A.GetKL();
    int ku = A.GetKU();
    for (int i = 0; i < A.GetM(); i++)
      for (int k = max(-i, -kl); k <= min(ku, this->n_-1-i); k++)
        {
          int j = i + k;
          AddInteraction(i, j, alpha*A(i, j));
        }
  }
  
  
  //! performs matrix-vector product y = y + alpha A x
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0, class T1>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::MltAdd(const T0& alpha, const SeldonTranspose& trans,
	   const Vector<T1>& x, Vector<T1>& y) const
  {
    int d = kl_ + ku_; T1 val;
    if (trans.NoTrans())
      {
        for (int j = 0; j < GetN(); j++)
          {
            int kmin = max(kl_, d - j), row;
            int kmax = min(kl_ + d, this->m_ -1 + d - j);
            val = alpha*x(j);
            for (int k = kmin; k <= kmax; k++)
              {
                row = j + k - d;
                y(row) = y(row) + data_(k, j)*val;
              }
          }
      }
    else if (trans.Trans())
      {
        for (int j = 0; j < GetN(); j++)
          {
            int kmin = max(kl_, d - j), row;
            int kmax = min(kl_ + d, this->m_ -1 + d - j);
            val = T1(0);
            for (int k = kmin; k <= kmax; k++)
              {
                row = j + k - d;
                val += data_(k, j)*x(row);
              }
            
            y(j) += alpha*val;
          }
      }
    else
      {
        for (int j = 0; j < GetN(); j++)
          {
            int kmin = max(kl_, d - j), row;
            int kmax = min(kl_ + d, this->m_ -1 + d - j);
            val = T1(0);
            for (int k = kmin; k <= kmax; k++)
              {
                row = j + k - d;
                val += conjugate(data_(k, j))*x(row);
              }
            
            y(j) += alpha*val;
          }
      }
  }
  
  
  //! solves A x = b, assuming that Factorize has been previously called
  template <class T, class Prop, class Storage, class Allocator>
  template<class T1>
  void Matrix_Band<T, Prop, Storage, Allocator>::Solve(Vector<T1>& x) const
  {
    int d = kl_ + ku_;
    // resolution of L y = x
    for (int j = 0; j < this->n_; j++)
      for (int p = 1; p <= min(this->m_ - 1 - j, kl_); p++)
        x(j+p) -= data_(d+p, j)*x(j);
    
    // resolution of U x = y
    for (int j = this->n_-1; j >= 0; j--)
      {        
        x(j) *= data_(d, j);
        for (int p = 1; p <= min(j, ku_); p++)
          x(j-p) -= data_(d-p, j)*x(j);
      }
  }


  //! solves A x = b, assuming that Factorize has been previously called
  //! (with pivoting)
  template <class T, class Prop, class Storage, class Allocator>
  template<class T1> void Matrix_Band<T, Prop, Storage, Allocator>
  ::Solve(const Vector<int>& ipivot, Vector<T1>& x) const
  {
    int d = kl_ + ku_;
    T1 tmp;
    // resolution of L y = x
    for (int j = 0; j < this->n_; j++)
      {
	// applying row interchange if needed
	if (ipivot(j) != j)
	  {
	    tmp = x(j);
	    x(j) = x(ipivot(j));
	    x(ipivot(j)) = tmp;
	  }
	
	for (int p = 1; p <= min(this->m_ - 1 - j, kl_); p++)
	  x(j+p) -= data_(d+p, j)*x(j);
      }
    
    // resolution of U x = y
    for (int j = this->n_-1; j >= 0; j--)
      {        
        x(j) *= data_(d, j);
        for (int p = 1; p <= min(j, kl_+ku_); p++)
          x(j-p) -= data_(d-p, j)*x(j);
      }    
  }
  
  
  //! writes the matrix in a file in binary format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Band::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();    
  }
  
  
  //! writes the matrix in binary format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->m_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->n_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->kl_)),
		     sizeof(int));
    FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&this->ku_)),
		     sizeof(int));
    
    data_.Write(FileStream, false);
  }
  
  
  //! writes the matrix in a file in text format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("TinyBandMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! writes the matrix in text format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Band<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {
    int d = kl_ + ku_;
    for (int j = 0; j < GetN(); j++)
      {
        int kmin = max(kl_, d - j);
	//int kmin = max(0, d-j);
        int kmax = min(kl_ + d, this->m_ -1 + d - j);
        for (int k = kmin; k <= kmax; k++)
          {
            int row = j + k - d;
            FileStream
	      << row + 1 << " " << j+1 << " " << data_(k, j) << '\n';
          }
      }            
  }
  
  
  /****************
   * Matrix_Arrow *
   ****************/
  

  //! clears the matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::Clear()
  {
    Matrix_Band<T, Prop, Storage, Allocator>::Clear();
    
    last_rows_.Clear();
    last_columns_.Clear();
    last_block_.Clear();
  }
  
  
  //! sets all non-zero entries to 0
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::Zero()
  {
    this->data_.Zero();
    last_rows_.Zero();
    last_columns_.Zero();
    last_block_.Zero();
  }
  
  
  //! changes the size of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::
  Reallocate(int m, int n, int kl, int ku,
             int nb_last_row, int nb_last_col)
  {
    Matrix_Band<T, Prop, Storage, Allocator>::
      Reallocate(m-nb_last_row, n-nb_last_col, kl, ku);
    
    last_rows_.Reallocate(nb_last_row, n-nb_last_col);
    last_columns_.Reallocate(m-nb_last_row, nb_last_col);
    last_block_.Reallocate(nb_last_row, nb_last_col);
  }


  //! performs the operation A(i, j) = A(i, j) + val
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::AddInteraction(int i, int j, const T& val)
  {
    if (i >= this->m_)
      {
        if (j >= this->n_)
          last_block_(i-this->m_, j-this->n_) += val;
        else
          last_rows_(i-this->m_, j) += val;
      }
    else
      {
        if (j >= this->n_)
          last_columns_(i, j-this->n_) += val;
        else
          Matrix_Band<T, Prop, Storage, Allocator>::
	    AddInteraction(i, j, val);
      }
  }
  
  
  //! adds severals values on a single row of the matrix  
  /*!
    \param[in] i row on which values are added
    \param[in] n number of values to add
    \param[in] num column numbers
    \param[in] val values
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::
  AddInteractionRow(int i, int n, const IVect& num, const Vector<T>& val,
		    bool sorted)
  {
    for (int j = 0; j < n; j++)
      AddInteraction(i, num(j), val(j));
  }
  
  
  //! clears row i, fills it with 0
  /*!
    \param[in] i row to clear
   */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::ClearRow(int i)
  {
    T zero; SetComplexZero(zero);
    if (i < this->m_)
      {
        for (int k = max(-i, -this->kl_);
	     k <= min(this->ku_, this->n_-1-i); k++)
          {
            int j = i+k;
            this->Get(i, j) = zero;
          }
        
        for (int j = 0; j < last_columns_.GetN(); j++)
          this->Get(i, this->n_+j) = zero;        
      }
    else
      {
        for (int j = 0; j < this->n_+this->last_columns_.GetN(); j++)
          this->Get(i, j) = zero;        
      }
  }
  
  
  //! returns entry (i, j) of the matrix
  template <class T, class Prop, class Storage, class Allocator>
  const T Matrix_Arrow<T, Prop, Storage, Allocator>
  ::operator()(int i, int j) const
  {
    if (i >= this->m_)
      {
        if (j >= this->n_)
          return last_block_(i-this->m_, j-this->n_);
        
        return last_rows_(i-this->m_, j);
      }

    if (j >= this->n_)
      return last_columns_(i, j-this->n_);
    
    return static_cast<const Matrix_Band<T, Prop, Storage, Allocator>& >
      (*this)(i, j);
  }
  
  
  //! multiplication by a scalar
  template <class T, class Prop, class Storage, class Allocator>
  Matrix<T, Prop, Storage, Allocator>& 
  Matrix_Arrow<T, Prop, Storage, Allocator>::operator *=(const T& alpha)
  {
    this->data_ *= alpha;
    last_rows_ *= alpha;
    last_columns_ *= alpha;
    last_block_ *= alpha;
    return static_cast<Matrix<T, Prop, Storage, Allocator>& >(*this);
  }

  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  T& Matrix_Arrow<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    if (i >= this->m_)
      {
        if (j >= this->n_)
          return last_block_(i-this->m_, j-this->n_);
        
        return last_rows_(i-this->m_, j);
      }

    if (j >= this->n_)
      return last_columns_(i, j-this->n_);
    
    return Matrix_Band<T, Prop, Storage, Allocator>::Get(i, j);
  }
  
  
  //! returns a reference to A(i, j) 
  template <class T, class Prop, class Storage, class Allocator>
  const T& Matrix_Arrow<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    if (i >= this->m_)
      {
        if (j >= this->n_)
          return last_block_(i-this->m_, j-this->n_);
        
        return last_rows_(i-this->m_, j);
      }

    if (j >= this->n_)
      return last_columns_(i, j-this->n_);
    
    return Matrix_Band<T, Prop, Storage, Allocator>::Get(i, j);
  }


  //! sets the matrix to the identity matrix
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::SetIdentity()
  {
    Matrix_Band<T, Prop, Storage, Allocator>::SetIdentity();
    
    T one, zero;
    SetComplexZero(zero);
    SetComplexOne(one);
    last_rows_.Fill(zero);
    last_columns_.Fill(zero);
    last_block_.Fill(zero);
    for (int i = min(this->m_, this->n_); i < this->GetM(); i++)
      this->Get(i, i) = one;
  }
    
  
  //! sets all non-zero entries to a same value
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::Fill(const T0& x)
  {
    T x_;
    SetComplexReal(x, x_);
    this->data_.Fill(x_);
    last_rows_.Fill(x_);
    last_columns_.Fill(x_);
    last_block_.Fill(x_);
  }
    
  
  //! sets all non-zero entries to random values
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::FillRand()
  {
    this->data_.FillRand();
    last_rows_.FillRand();
    last_columns_.FillRand();
    last_block_.FillRand();
  }
  
  
  //! performs LU factorisation without pivoting
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::Factorize()
  {
    // position of the main diagonal
    int d = this->kl_ + this->ku_;
    
    // treating first columns
    T pivot, diag;
    int nb_last_row = last_rows_.GetM();
    int nb_last_col = last_columns_.GetN();
    T one; SetComplexOne(one);
    for (int j = 0; j < this->n_; j++)
      {
        // replacing diagonal element by its inverse
	if (j >= this->m_)
	  {
	    diag = one / this->last_rows_(j-this->m_, j);
	    this->last_rows_(j-this->m_, j) = diag;
	  }
	else
	  {
	    this->data_(d, j) = one / this->data_(d, j);
	    diag = this->data_(d, j);
	  }
	
	// Gaussian elimination of rows below row j (part in the band)
        for (int p = 1; p <= min(this->m_ - 1 - j, this->kl_); p++)
          {
            pivot = this->data_(d+p, j)*diag;
            this->data_(d+p, j) = pivot;
            // band
            for (int k = 1; k <= min(this->n_ - 1 - j, this->ku_); k++)
              this->data_(d+p-k, j+k) -= pivot*this->data_(d-k, j+k);
            
            // treating coefficients of the last columns
            for (int k = 0; k < nb_last_col; k++)
              last_columns_(j+p, k) -= pivot*last_columns_(j, k);
          }

        // treating coefficients of the last rows
        for (int p = max(0, j-this->m_+1); p < nb_last_row; p++)
          {
            pivot = last_rows_(p, j)*diag;
            last_rows_(p, j) = pivot;
            if (j >= this->m_)
	      {
		for (int k = j+1; k < this->n_; k++)
		  last_rows_(p, k) -= pivot*this->last_rows_(j-this->m_, k);

		// treating coefficients of the last block
		for (int k = 0; k < nb_last_col; k++)
		  last_block_(p, k) -= pivot*last_block_(j-this->m_, k);
	      }
	    else
	      {
		for (int k = 1; k <= min(this->n_ - 1 - j, this->ku_); k++)
		  last_rows_(p, j+k) -= pivot*this->data_(d-k, j+k);
		
		// treating coefficients of the last block
		for (int k = 0; k < nb_last_col; k++)
		  last_block_(p, k) -= pivot*last_columns_(j, k);
	      }
          }
      }
    
    // treating last columns
    for (int j = this->n_; j < this->GetM(); j++)
      {
        if (j >= this->m_)
	  {
	    diag = one / last_block_(j-this->m_, j-this->n_);
	    last_block_(j-this->m_, j-this->n_) = diag;
	  }
	else
	  {
	    diag = one / last_columns_(j, j-this->n_);
	    last_columns_(j, j-this->n_) = diag;
	    
	    // elimination in rows below but still located in last_columns
	    for (int i = j+1; i < this->m_; i++)
	      {
		pivot = last_columns_(i, j-this->n_)*diag;
		last_columns_(i, j-this->n_) = pivot;
		for (int k = j+1; k < this->GetM(); k++)
		  last_columns_(i, k-this->n_) 
		    -= pivot*last_columns_(j, k-this->n_);
	      }
	  }
	
	// now elimination of rows in last_block
        for (int i = max(j+1,this->m_); i < this->GetM(); i++)
          {
            pivot = last_block_(i-this->m_, j-this->n_)*diag;
            last_block_(i-this->m_, j-this->n_) = pivot;
            if (j >= this->m_)
	      for (int k = j+1; k < this->GetM(); k++)
		last_block_(i-this->m_, k-this->n_) 
		  -= pivot*last_block_(j-this->m_, k-this->n_);
	    else
	      for (int k = j+1; k < this->GetM(); k++)
		last_block_(i-this->m_, k-this->n_) 
		  -= pivot*last_columns_(j, k-this->n_);
          }
      }    
  }
    
  
  //! performs the operation *this = *this + alpha A
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::
  Add_(const T0& alpha, const Matrix<T, General, ArrowCol, Allocator>& A)
  {
    // adding main band of A
    int kl = A.GetKL();
    int ku = A.GetKU();
    int m = A.GetM() - A.GetNbLastRow();
    int n = A.GetN() - A.GetNbLastCol();
    for (int i = 0; i < m; i++)
      for (int k = max(-i, -kl); k <= min(ku, n-1-i); k++)
        {
          int j = i + k;
          AddInteraction(i, j, alpha*A(i, j));
        }
    
    // then last rows
    for (int i = 0; i < A.GetNbLastRow(); i++)
      for (int j = 0; j < A.GetN(); j++)
        AddInteraction(m+i, j, alpha*A(m+i, j));
    
    // and last columns
    for (int j = 0; j < A.GetNbLastCol(); j++)
      for (int i = 0; i < m; i++)
        AddInteraction(i, n+j, alpha*A(i, n+j));
  }
  

  //! performs matrix-vector product y = y + alpha A x  
  template <class T, class Prop, class Storage, class Allocator>
  template<class T0, class T1>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::
  MltAdd(const T0& alpha, const SeldonTranspose& trans,
	 const Vector<T1>& x, Vector<T1>& y) const
  {
    // banded part
    Matrix_Band<T, Prop, Storage, Allocator>::MltAdd(alpha, trans, x, y);
    
    T1 val, zero; SetComplexZero(zero);
    if (trans.NoTrans())
      {
        // last rows
        for (int i = 0; i < last_rows_.GetM(); i++)
          {
            val = zero;
            for (int j = 0; j < this->n_; j++)
              val += last_rows_(i, j)*x(j);

            y(this->m_ + i) += alpha*val;
          }
        
        // last columns
        for (int j = 0; j < last_columns_.GetN(); j++)
          {
            val = alpha*x(this->n_+j);
            for (int i = 0; i < this->m_; i++)
              y(i) += last_columns_(i, j)*val;
          }
        
        // last block
        for (int i = 0; i < last_block_.GetM(); i++)
          {
            val = zero;
            for (int j = 0; j < last_block_.GetN(); j++)
              val += last_block_(i, j)*x(this->n_+j);
            
            y(this->m_+i) += alpha*val;
          }
      }
    else if (trans.Trans())
      {
        // last rows
        for (int i = 0; i < last_rows_.GetM(); i++)
          {
            val = alpha*x(this->m_ + i);
            for (int j = 0; j < this->n_; j++)
              y(j) += last_rows_(i, j)*val;
          }
        
        // last columns
        for (int j = 0; j < last_columns_.GetN(); j++)
          {
            val = zero;
            for (int i = 0; i < this->m_; i++)
              val += last_columns_(i, j)*x(i);
        
            y(this->n_+j) += alpha*val;
          }
        
        // last block
        for (int i = 0; i < last_block_.GetM(); i++)
          {
            val = alpha*x(this->m_+i);
            for (int j = 0; j < last_block_.GetN(); j++)
              y(this->n_+j) += last_block_(i, j)*val;
          }
      }
    else
      {
        // last rows
        for (int i = 0; i < last_rows_.GetM(); i++)
          {
            val = alpha*x(this->m_ + i);
            for (int j = 0; j < this->n_; j++)
              y(j) += conjugate(last_rows_(i, j))*val;
          }
        
        // last columns
        for (int j = 0; j < last_columns_.GetN(); j++)
          {
            val = zero;
            for (int i = 0; i < this->m_; i++)
              val += conjugate(last_columns_(i, j))*x(i);
        
            y(this->n_+j) += alpha*val;
          }
        
        // last block
        for (int i = 0; i < last_block_.GetM(); i++)
          {
            val = alpha*x(this->m_+i);
            for (int j = 0; j < last_block_.GetN(); j++)
              y(this->n_+j) += conjugate(last_block_(i, j))*val;
          }
      }
  }
  
  
  //! solves A x = b, assuming that Factorize has been previously called
  template <class T, class Prop, class Storage, class Allocator>
  template<class T1>
  void Matrix_Arrow<T, Prop, Storage, Allocator>::Solve(Vector<T1>& x) const
  {
    int d = this->kl_ + this->ku_;
    // resolution of L y = x
    
    // part due to band
    for (int j = 0; j < min(this->m_, this->n_); j++)
      for (int p = 1; p <= min(this->m_ - 1 - j, this->kl_); p++)
        x(j+p) -= this->data_(d+p, j)*x(j);
    
    // part between the band and the last block
    if (this->m_ > this->n_)
      for (int j = this->n_; j < this->m_; j++)
	for (int k = j+1; k < this->m_; k++)
	  x(k) -= last_columns_(k, j-this->n_)*x(j);
    else
      for (int j = this->m_; j < this->n_; j++)
	for (int k = 0; k < j; k++)
	  x(j) -= last_rows_(j-this->m_, k)*x(k);
    
    // part due to the last block
    for (int i = max(this->m_, this->n_); i < this->GetM(); i++)
      {
        for (int j = 0; j < this->n_; j++)
          x(i) -= last_rows_(i-this->m_, j)*x(j);
        
        for (int j = this->n_; j < i; j++)
          x(i) -= last_block_(i-this->m_, j-this->n_)*x(j);
      }
    
    // resolution of U x = y
    
    // part due to the last block
    int N = this->GetM();
    for (int i = N-1; i >= max(this->m_, this->n_); i--)
      {
	x(i) *= last_block_(i-this->m_, i-this->n_);
        for (int j = this->m_; j < i; j++)
          x(j) -= last_block_(j-this->m_, i-this->n_)*x(i);
	
	for (int j = 0; j < this->m_; j++)
	  x(j) -= last_columns_(j, i-this->n_)*x(i);
      }
    
    // part between band and the last block
    if (this->m_ > this->n_)
      for (int i = this->m_-1; i >= this->n_; i--)
	{
	  x(i) *= last_columns_(i, i-this->n_);
	  for (int j = 0; j < i; j++)
	    x(j) -= last_columns_(j, i-this->n_)*x(i);
	}
    else
      for (int i = this->n_-1; i >= this->m_; i--)
	{
	  for (int j = i+1; j < this->n_; j++)
	    x(i) -= last_rows_(i-this->m_, j)*x(j);
	  
	  x(i) *= last_rows_(i-this->m_, i);
	}
      
    // part due to the band
    for (int j = this->n_-1; j >= 0; j--)
      {        
        if (j < min(this->m_, this->n_))
	  x(j) *= this->data_(d, j);
	
        for (int p = 1; p <= min(j, this->ku_); p++)
          if (j-p < this->m_)
	    x(j-p) -= this->data_(d-p, j)*x(j);
      }
  }
    

  //! writes the matrix in a file in binary format  
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::Write(string FileName) const
  {
    ofstream FileStream;
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("Matrix_Arrow::Write(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->Write(FileStream);

    FileStream.close();    
  }
  
  
  //! writes the matrix in binary format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::Write(ostream& FileStream) const
  {
    Matrix_Band<T, Prop, Storage, Allocator>::Write(FileStream);
    
    last_rows_.Write(FileStream);
    last_columns_.Write(FileStream);
    last_block_.Write(FileStream);
  }
  
  
  //! writes the matrix in a file in text format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::WriteText(string FileName) const
  {
    ofstream FileStream;
    FileStream.precision(cout.precision());
    FileStream.flags(cout.flags());
    FileStream.open(FileName.c_str());

#ifdef SELDON_CHECK_IO
    // Checks if the file was opened.
    if (!FileStream.is_open())
      throw IOError("TinyBandMatrix::WriteText(string FileName)",
		    string("Unable to open file \"") + FileName + "\".");
#endif

    this->WriteText(FileStream);

    FileStream.close();
  }
  
  
  //! writes the matrix in text format
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_Arrow<T, Prop, Storage, Allocator>
  ::WriteText(ostream& FileStream) const
  {
    int d = this->kl_ + this->ku_;
    for (int j = 0; j < this->n_; j++)
      {
        int kmin = max(this->kl_, d - j);
        int kmax = min(this->kl_ + d, this->m_ -1 + d - j);
        for (int k = kmin; k <= kmax; k++)
          {
            int row = j + k - d;
            FileStream
	      << row + 1 << " " << j+1 << " " << this->data_(k, j) << '\n';
          }
        
        for (int k = 0; k < last_rows_.GetM(); k++)
          {
            int row = this->m_ + k;
            FileStream
	      << row + 1 << " " << j+1 << " " << last_rows_(k, j) << '\n';
          }
      }    
    
    for (int j = 0; j < last_columns_.GetN(); j++)
      {
        for (int i = 0; i < this->m_; i++)
          FileStream << i+1 << " " << this->n_+j+1 << " "
		     << last_columns_(i, j) << '\n';
        
        for (int i = 0; i < last_block_.GetM(); i++)
          FileStream << this->m_+i+1 << " " << this->n_+j+1
		     << " " << last_block_(i, j) << '\n';
      }
  }
  
  
  /*************
   * Functions *
   *************/


  //! LU factorisation
  template<class T, class Allocator>
  void GetLU(Matrix<T, General, BandedCol, Allocator>& A,
             Matrix<T, General, BandedCol, Allocator>& mat_lu,
	     bool keep_matrix)
  {
    mat_lu = A;
    if (!keep_matrix)
      A.Clear();
    
    mat_lu.Factorize();
  }

  
  //! LU factorisation
  template<class T, class Allocator>
  void GetLU(Matrix<T, General, BandedCol, Allocator>& A)
  {
    A.Factorize();
  }
  

  //! LU factorisation with row interchanges
  template<class Allocator>
  void GetLU(Matrix<double, General, BandedCol, Allocator>& A,
             Vector<int>& ipivot, LapackInfo& info)
  {
#ifdef SELDON_WITH_LAPACK
    int m = A.GetM();
    int n = A.GetN();
#ifdef SELDON_CHECK_BOUNDS
    if ((m <= 0)||(n <= 0))
      throw WrongDim("GetLU", "Provide a non-empty matrix");
#endif
    
    ipivot.Reallocate(m);
    int kl = A.GetKL();
    int ku = A.GetKU();
    int lda = 2*kl + ku + 1;
    dgbtrf_(&m, &n, &kl, &ku, A.GetData(), &lda,
            ipivot.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetLU",
			"An error occured during the factorization.");
#endif

#else
    
    A.Factorize(ipivot);
    
#endif
  }
  
  
  //! resolution of A x = b, assuming that GetLU has been called
  template<class Allocator>
  void SolveLU(const Matrix<double, General, BandedCol, Allocator>& A,
               const Vector<int>& ipivot, Vector<double>& b,
               LapackInfo& info)
  {
#ifdef SELDON_WITH_LAPACK
    int n = A.GetM();
    int kl = A.GetKL();
    int ku = A.GetKU();
    int nrhs = 1;
    int lda = 2*kl + ku + 1;
    int ldb = n;
    char trans('N');
    dgbtrs_(&trans, &n, &kl, &ku, &nrhs, A.GetData(), &lda,
            ipivot.GetData(), b.GetData(), &ldb, &info.GetInfoRef());
#else

    A.Solve(ipivot, b);
#endif
    
  }  
  
  
  //! resolution of A x = b, assuming that GetLU has been called
  template<class Allocator>
  void SolveLU(const Matrix<double, General, BandedCol, Allocator>& A,
               const Vector<int>& ipivot, Vector<complex<double> >& b,
               LapackInfo& info)
  {
#ifdef SELDON_WITH_LAPACK
    int n = A.GetM();
    int kl = A.GetKL();
    int ku = A.GetKU();
    int nrhs = 2;
    int lda = 2*kl + ku + 1;
    int ldb = n;
    char trans('N');

    Matrix<double, General, ColMajor> brhs(n, 2);
    for (int i = 0; i < n; i++)
      {
        brhs(i, 0) = real(b(i));
        brhs(i, 1) = imag(b(i));
      }
    
    dgbtrs_(&trans, &n, &kl, &ku, &nrhs, A.GetData(), &lda,
            ipivot.GetData(), brhs.GetData(), &ldb, &info.GetInfoRef());
    
    for (int i = 0; i < n; i++)
      b(i) = complex<double>(brhs(i, 0), brhs(i, 1));

#else
    
    A.Solve(ipivot, b);
#endif

  }  


  //! LU factorisation with row interchanges
  template<class Allocator>
  void GetLU(Matrix<complex<double>, General, BandedCol, Allocator>& A,
             Vector<int>& ipivot, LapackInfo& info)
  {
#ifdef SELDON_WITH_LAPACK
    int m = A.GetM();
    int n = A.GetN();
#ifdef SELDON_CHECK_BOUNDS
    if ((m <= 0)||(n <= 0))
      throw WrongDim("GetLU", "Provide a non-empty matrix");
#endif

    ipivot.Reallocate(m);    
    int kl = A.GetKL();
    int ku = A.GetKU();
    int lda = 2*kl + ku + 1;
    zgbtrf_(&m, &n, &kl, &ku, A.GetData(), &lda,
            ipivot.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetLU",
			"An error occured during the factorization.");
#endif

#else
    
    A.Factorize(ipivot);    
#endif

  }
  
  
  //! resolution of A x = b, assuming that GetLU has been called
  template<class Allocator>
  void SolveLU(const Matrix<complex<double>,
	       General, BandedCol, Allocator>& A,
               const Vector<int>& ipivot, Vector<complex<double> >& b,
               LapackInfo& info)
  {
#ifdef SELDON_WITH_LAPACK
    int n = A.GetM();
    int kl = A.GetKL();
    int ku = A.GetKU();
    int nrhs = 1;
    int lda = 2*kl + ku + 1;
    int ldb = n;
    char trans('N');
    zgbtrs_(&trans, &n, &kl, &ku, &nrhs, A.GetData(), &lda,
            ipivot.GetData(), b.GetDataVoid(), &ldb, &info.GetInfoRef());
#else
    
    A.Solve(ipivot, b);
#endif

  }  
  
  
  //! factorizes matrix with pivoting
  template<class T, class Allocator>
  void GetLU(Matrix<T, General, BandedCol, Allocator>& A,
             Vector<int>& ipivot, LapackInfo& info)
  {
    A.Factorize(ipivot);
  }


  //! Solves A x = b with pivoting
  template<class T, class Allocator>
  void SolveLU(const Matrix<T, General, BandedCol, Allocator>& A,
	       const Vector<int>& ipivot, Vector<T>& b, LapackInfo& info)
  {
    A.Solve(ipivot, b);
  }


  //! Solves A x = b with pivoting
  template<class T, class Allocator>
  void SolveLU(const Matrix<T, General, BandedCol, Allocator>& A,
	       const Vector<int>& ipivot, Vector<complex<T> >& b, LapackInfo& info)
  {
    A.Solve(ipivot, b);
  }


  //! Scaling of a matrix
  /*!
    A is replaced by Drow A Dcol
    where Drow and Dcol are diagonal matrices
    and stored as dense vectors
  */
  template<class T, class Prop, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(Matrix<T, Prop, BandedCol, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol)
  {
    int kl = A.GetKL(), ku = A.GetKU(), n = A.GetM();
    for (int i = 0; i < A.GetM(); i++)
      for (int j = max(0, i-kl); j < min(n, i+ku+1); j++)
        A.Get(i, j) *= Drow(i)*Dcol(j);    
  }
  
}

#define SELDON_FILE_BAND_MATRIX_CXX
#endif
