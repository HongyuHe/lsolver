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


#ifndef SELDON_FILE_COMMON_INLINE_CXX


namespace std
{
  template<class T>
  inline T conjugate(const T& x)
  {
    return x;
  }

  template<class T>
  inline complex<T> conjugate(const complex<T>& x)
  {
    return conj(x);
  }

  template<class T>
  inline T realpart(const T& x)
  {
    return x;
  }

  template<class T>
  inline T realpart(const complex<T>& x)
  {
    return real(x);
  }
}

namespace Seldon
{

  template <class T>
  inline void PrintArray(T* v, int lgth)
  {
    for (int k = 0; k < lgth - 1; k++)
      std::cout << v[k] << " | ";
    std::cout << v[lgth - 1] << std::endl;
  }


  //! Default constructor.
  inline Str::Str()
  {
  }


  //! Copy constructor.
  /*!
    \param[in] s 'Str' instance to be copied.
  */
  inline Str::Str(const Str& s)
  {
    string s_input = s;
    output_ << s_input;
  }


  //! Conversion to string.
  inline Str::operator std::string() const
  {
    return output_.str();
  }


  //! Adds an element to the string.
  /*!
    \param[in] input element added at the end of the string.
  */
  template <class T>
  inline Str& Str::operator << (const T& input)
  {
    output_ << input;
    return *this;
  }


  //! Adds an element to an instance of 'Str'.
  /*!
    \param[in] s 'Str' instance.
    \param[in] input element added at the end of the string.
  */
  template <class T>
  inline Str operator + (const Str& s, const T& input)
  {
    string s_input = s;
    Str output;
    output << s_input << input;
    return output;
  }


  //! Converts a 'str' instance to an 'ostream' instance.
  inline ostream& operator << (ostream& out, Str& in)
  {
    string output = in;
    out << output;
    return out;
  }


  //! Converts a 'str' instance to an 'ostream' instance.
  inline ostream& operator << (ostream& out, Str in)
  {
    string output = in;
    out << output;
    return out;
  }


  //! Converts most types to string.
  /*!
    \param input variable to be converted.
    \return A string containing 'input'.
  */
  template<typename T>
  inline std::string to_str(const T& input)
  {
    std::ostringstream output;
    output.precision(cout.precision());
    output << input;
    return output.str();
  }


  //! Converts string to most types, specially numbers.
  /*!
    \param[in] s string to be converted.
    \param[out] num \a s converted to 'T'.
  */
  template <class T>
  inline void to_num(std::string s, T& num)
  {
    std::istringstream str(s);
    str >> num;
  }


  //! Converts string to most types, specially numbers.
  /*!
    \param[in] s string to be converted.
    \return \a s converted to 'T'.
  */
  template <class T>
  inline T to_num(std::string s)
  {
    T num;
    std::istringstream str(s);
    str >> num;
    return num;
  }


  //! Sets a number to zero.
  /*!
    \param[in,out] number number to be set to zero.
  */
  template <class T>
  inline void SetComplexZero(T& number)
  {
    number = T(0);
  }


  //! Sets a complex number to zero.
  /*!
    \param[in,out] number complex number to be set to zero.
  */
  template <class T>
  inline void SetComplexZero(std::complex<T>& number)
  {
    number = std::complex<T>(T(0), T(0));
  }


  //! Sets a number to one.
  /*!
    \param[in,out] number number to be set to one.
  */
  template <class T>
  inline void SetComplexOne(T& number)
  {
    number = T(1);
  }


  //! Sets a complex number to (1, 0).
  /*!
    \param[in,out] number complex number to be set to (1, 0).
  */
  template <class T>
  inline void SetComplexOne(std::complex<T>& number)
  {
    number = complex<T>(T(1), T(0));
  }


  //! Sets a real number to n.
  /*!
    \param[in,out] number real umber to be set to n.
  */
  template <class T>
  inline void SetComplexReal(int n, T& number)
  {
    number = T(n);
  }


  //! Sets a complex number to (n, 0).
  /*!
    \param[in,out] number complex number to be set to (n, 0).
  */
  template <class T>
  inline void SetComplexReal(int n, std::complex<T>& number)
  {
    number = std::complex<T>(n, 0);
  }


  //! Sets a complex number to (n, 0).
  /*!
    \param[in,out] number complex number to be set to (n, 0).
  */
  template <class T>
  inline void SetComplexReal(bool n, std::complex<T>& number)
  {
    number = std::complex<T>(n, 0);
  }


  //! Sets a complex number to (n, 0).
  inline void SetComplexReal(int x, std::complex<int>& number)
  {
    number = std::complex<int>(x, 0);
  }


  //! Sets a complex number to (x, 0).
  /*!
    \param[in,out] number complex number to be set to (x, 0).
  */
  template <class T>
  inline void SetComplexReal(const T& x, std::complex<T>& number)
  {
    number = std::complex<T>(x, 0);
  }


  //! Sets a complex number to (x, 0).
  /*!
    \param[in,out] number complex number to be set to (x, 0).
  */
  template <class T>
  inline void SetComplexReal(const std::complex<T>& x, T& number)
  {
    number = std::realpart(x);
  }


  //! Sets a complex number to x.
  /*!
    \param[in,out] number complex number to be set to x.
  */
  template <class T0, class T1>
  inline void SetComplexReal(const T0& x, T1& number)
  {
    number = x;
  }


  //! Returns true for a complex number
  template <class T>
  inline bool IsComplexNumber(const T& number)
  {
    return false;
  }


  //! Returns true for a complex number
  template <class T>
  inline bool IsComplexNumber(const std::complex<T>& number)
  {
    return true;
  }


  //! returns absolute value of val
  template<class T>
  inline T ComplexAbs(const T& val)
  {
    return abs(val);
  }


  //! returns modulus of val
  template<class T>
  inline T ComplexAbs(const std::complex<T>& val)
  {
#if defined(SELDON_WITH_BLAS) && !defined(SELDON_WITH_LAPACK)
    // we choose Blas convention
    return abs(real(val)) + abs(imag(val));
#else
    // otherwise usual modulus
    return abs(val);
#endif
  }


  //! returns the square modulus of z
  template<class T>
  inline T absSquare(const std::complex<T>& z)
  {
    // more optimal than real(z * conj(z))
    return real(z)*real(z) + imag(z)*imag(z);
  }


  //! returns the square modulus of z
  template<class T>
  inline T absSquare(const T& z)
  {
    return z*z;
  }

}  // namespace Seldon.

#define SELDON_FILE_COMMON_INLINE_CXX
#endif
