// Copyright (C) 2010 Vivien Mallet
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


#ifndef SELDON_WITH_BLAS
#define SELDON_WITH_BLAS
#endif

#include "SeldonHeader.hxx"

#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/MatrixFlag.cxx"
#include "share/CommonInline.cxx"
#include "share/StorageInline.cxx"
#include "share/Errors.cxx"
#endif

#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;


namespace Seldon
{

  SELDON_EXTERN template std::string to_str<int>(const int&);
  SELDON_EXTERN template std::string to_str<long>(const long&);
  SELDON_EXTERN template std::string to_str<float>(const float&);
  SELDON_EXTERN template std::string to_str<double>(const double&);
  SELDON_EXTERN template std::string to_str<complex<float> >(const complex<float>&);
  SELDON_EXTERN template std::string to_str<complex<double> >(const complex<double>&);

  SELDON_EXTERN template void to_num(std::string, int&);
  SELDON_EXTERN template void to_num(std::string, float&);
  SELDON_EXTERN template void to_num(std::string, double&);
  SELDON_EXTERN template void to_num(std::string, complex<float>&);
  SELDON_EXTERN template void to_num(std::string, complex<double>&);

  SELDON_EXTERN template int to_num<int>(std::string);
  SELDON_EXTERN template float to_num<float>(std::string);
  SELDON_EXTERN template double to_num<double>(std::string);
  SELDON_EXTERN template complex<float> to_num<complex<float> >(std::string);
  SELDON_EXTERN template complex<double> to_num<complex<double> >(std::string);

  SELDON_EXTERN template void SetComplexZero(@scalar&);
  SELDON_EXTERN template void SetComplexOne(@scalar&);

}
