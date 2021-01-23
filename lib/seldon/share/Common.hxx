// Copyright (C) 2001-2012 Vivien Mallet
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


#ifndef SELDON_FILE_COMMON_HXX

#include <complex>
#include <iostream>
#include <typeinfo>

#ifdef SELDON_WITH_HDF5
#include <hdf5.h>
#endif

namespace std
{
  template<class T>
  T conjugate(const T& x);

  template<class T>
  complex<T> conjugate(const complex<T>& x);

  template<class T>
  T realpart(const T& x);

  template<class T>
  T realpart(const complex<T>& x);  
}
      
namespace Seldon
{
  using namespace std;
 
  template <class T>
  void PrintArray(T* v, int lgth);
  
  //! This class helps formatting C++ strings on the fly.
  /*!
    It should may be used like that:
    string output = Str() + "There are " + 3 + " laws of robotics.";
  */
  class Str
  {
  private:
    //! Buffer.
    std::ostringstream output_;

  public:
    Str();
    Str(const Str& s);
    operator std::string() const;
    template <class T>
    Str& operator << (const T& input);
  };

  template <class T>
  Str operator + (const Str&, const T& input);

#ifndef SWIG
  ostream& operator << (ostream& out, Str& in);
  ostream& operator << (ostream& out, Str in);
#endif


  template<typename T>
  std::string to_str(const T& input);

  template <class T>
  void to_num(std::string s, T& num);

  template <class T>
  T to_num(std::string s);

  //! workaround class to retrieve double type from complex<double>
  /*!
    when T is a template parameter (double or complex<double>)
    if you need to know the associated real type (double),
    you can type :
    typedef typename ClassComplexType<T>::Treal real;
    In the same way, you can get complex<double> type :
    typedef typename ClassComplexType<T>::Tcplx complexe;
   */
  template<class T>
  class ClassComplexType
  {
    public : 
    typedef T Treal;
    typedef std::complex<T> Tcplx;
  };

  //! workaround class to retrieve double type from complex<double>
  template<class T>
  class ClassComplexType< complex<T> >
  {
  public : 
    typedef T Treal;
    typedef std::complex<T> Tcplx;
  };
  
  template <class T>
  void SetComplexZero(T& number);

  template <class T>
  void SetComplexZero(std::complex<T>& number);

  template <class T>
  void SetComplexOne(T& number);

  template <class T>
  void SetComplexOne(std::complex<T>& number);

  template <class T>
  void SetComplexReal(int n, std::complex<T>& number);

  template <class T>
  void SetComplexReal(int n, T& number);

  template <class T>
  void SetComplexReal(bool n, std::complex<T>& number);

  void SetComplexReal(int x, std::complex<int>& number);

  template <class T>
  void SetComplexReal(const T& x, std::complex<T>& number);

  template <class T>
  void SetComplexReal(const std::complex<T>& x, T& number);
  
  template <class T0, class T1>
  void SetComplexReal(const T0& x, T1& number);
  
  template<class T>
  T ComplexAbs(const T& val);

  template<class T>
  T ComplexAbs(const std::complex<T>& val);

  template<class T>
  T absSquare(const T& x);

  template<class T>
  T absSquare(const std::complex<T>& x);
  
  string GetExtension(const string& nom);
  string GetBaseString(const string& nom);
  
#ifdef SELDON_WITH_HDF5
  template <class T>
  hid_t GetH5Type(T& input);
#endif

  void InitSeldon(int, char**);
  int FinalizeSeldon();
  
}  // namespace Seldon.

#include "CommonInline.cxx"
#define SELDON_FILE_COMMON_HXX
#endif
