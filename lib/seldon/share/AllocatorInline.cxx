// Copyright (C) 2001-2010 Vivien Mallet, INRIA
// Author(s): Vivien Mallet, Marc Fragu
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


#ifndef SELDON_FILE_ALLOCATOR_INLINE_CXX

#include "Allocator.hxx"

namespace Seldon
{


  /////////////////
  // MALLOCALLOC //
  /////////////////


  template <class T>
  inline typename MallocAlloc<T>::pointer
  MallocAlloc<T>::allocate(size_t num, void* h)
  {
    return static_cast<pointer>( malloc(num * sizeof(T)) );
  }

  template <class T>
  inline void MallocAlloc<T>::deallocate(pointer data, size_t num, void* h)
  {
    free(data);
  }

  template <class T>
  inline void* MallocAlloc<T>::reallocate(pointer data, size_t num, void* h)
  {
    return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
  }

  template <class T>
  inline void MallocAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  MallocAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }


  /////////////////
  // CALLOCALLOC //
  /////////////////


  template <class T>
  inline typename CallocAlloc<T>::pointer
  CallocAlloc<T>::allocate(size_t num, void* h)
  {
    return static_cast<pointer>( calloc(num, sizeof(T)) );
  }

  template <class T>
  inline void CallocAlloc<T>::deallocate(pointer data, size_t num, void* h)
  {
    free(data);
  }

  template <class T>
  inline void* CallocAlloc<T>::reallocate(pointer data, size_t num, void* h)
  {
    return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
  }

  template <class T>
  inline void CallocAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  CallocAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }


  //////////////
  // NEWALLOC //
  //////////////


  template <class T>
  inline typename NewAlloc<T>::pointer
  NewAlloc<T>::allocate(size_t num, void* h)
  {
    return static_cast<pointer>(new T[num]);
  }

  template <class T>
  inline void NewAlloc<T>::deallocate(pointer data, size_t num, void* h)
  {
    delete [] data;
  }

  template <class T>
  inline void* NewAlloc<T>::reallocate(pointer data, size_t num, void* h)
  {
    if (data != NULL)
      delete [] data;
    
    return (new T[num]);
  }

  template <class T>
  inline void NewAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  NewAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    for (size_t i = 0; i < num; i++)
      datat[i] = datas[i];
  }


  //////////////
  // NANALLOC //
  //////////////


  template <class T>
  inline void NaNAlloc<T>::deallocate(pointer data, size_t num, void* h)
  {
    free(data);
  }

  template <class T>
  inline void NaNAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  NaNAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }


} // namespace Seldon.

#define SELDON_FILE_ALLOCATOR_INLINE_CXX
#endif
