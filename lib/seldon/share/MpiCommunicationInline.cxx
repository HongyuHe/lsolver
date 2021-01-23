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

#ifndef SELDON_FILE_MPI_COMMUNICATION_INLINE_CXX

namespace Seldon
{

  inline MPI_Datatype GetMpiDataType(const Vector<bool>&)
  {
    return MPI_BYTE;
  }

  inline MPI_Datatype GetMpiDataType(const Vector<int>&)
  {
    return MPI_INTEGER;
  }

  inline MPI_Datatype GetMpiDataType(const Vector<long>&)
  {
    return MPI_LONG;
  }
  
  inline MPI_Datatype GetMpiDataType(const Vector<float>&)
  {
    return MPI_FLOAT;
  }
  
  inline MPI_Datatype GetMpiDataType(const Vector<complex<float> >&)
  {
    return MPI_FLOAT;
  }
  
  inline MPI_Datatype GetMpiDataType(const Vector<double>&)
  {
    return MPI_DOUBLE;
  }
  
  inline MPI_Datatype GetMpiDataType(const Vector<complex<double> >&)
  {
    return MPI_DOUBLE;
  }

  inline MPI_Datatype GetMpiDataType(const Vector<long double>&)
  {
    return MPI_LONG_DOUBLE;
  }
  
  inline MPI_Datatype GetMpiDataType(const Vector<complex<long double> >&)
  {
    return MPI_LONG_DOUBLE;
  }

  inline MPI_Datatype GetMpiDataType(const int&)
  {
    return MPI_INTEGER;
  }

  inline MPI_Datatype GetMpiDataType(const long&)
  {
    return MPI_LONG;
  }

  inline MPI_Datatype GetMpiDataType(const bool&)
  {
    return MPI_BYTE;
  }
  
  inline MPI_Datatype GetMpiDataType(const float&)
  {
    return MPI_FLOAT;
  }
  
  inline MPI_Datatype GetMpiDataType(const complex<float>&)
  {
    return MPI_FLOAT;
  }
  
  inline MPI_Datatype GetMpiDataType(const double&)
  {
    return MPI_DOUBLE;
  }
  
  inline MPI_Datatype GetMpiDataType(const complex<double>&)
  {
    return MPI_DOUBLE;
  }

  inline MPI_Datatype GetMpiDataType(const long double&)
  {
    return MPI_LONG_DOUBLE;
  }
  
  inline MPI_Datatype GetMpiDataType(const complex<long double>&)
  {
    return MPI_LONG_DOUBLE;
  }
  
  template<class T>
  inline int GetRatioMpiDataType(const Vector<T>&)
  {
    return 1;
  }
  
  template<class T>
  inline int GetRatioMpiDataType(const Vector<complex<T> >&)
  {
    return 2;
  }

  template<class T>
  inline int GetRatioMpiDataType(const T&)
  {
    return 1;
  }
  
  template<class T>
  inline int GetRatioMpiDataType(complex<T>&)
  {
    return 2;
  }
  
}

#define SELDON_FILE_MPI_COMMUNICATION_INLINE_CXX
#endif
