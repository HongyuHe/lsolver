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


#ifndef SELDON_FILE_ALLOCATOR_CXX

#include "Allocator.hxx"

namespace Seldon
{


  //////////////////
  // MALLOCOBJECT //
  //////////////////


  template <class T>
  typename MallocObject<T>::pointer
  MallocObject<T>::allocate(size_t num, void* h)
  {
    // The cast from char* to T* may lead to a memory shift (because of
    // alignment issues) under MS Windows. It requires that one allocates more
    // memory than necessary for memory_block.
    
    void* memory_block = malloc(sizeof(size_t) + sizeof(char*) +
                                (num + 2) * sizeof(T));
    memcpy(memory_block, &num, sizeof(size_t));
    char* data = static_cast<char*>(memory_block)
      + sizeof(size_t) + sizeof(char*) + sizeof(T);

    // The memory shift can occur here.
    pointer data_P = reinterpret_cast<pointer>(new(data) T[num]);
    size_t taille = sizeof(char*);
    memcpy(reinterpret_cast<char *>(data_P) - sizeof(char*),
           &memory_block, taille);
           
    return data_P;
  }


  template <class T>
  void MallocObject<T>::deallocate(pointer data, size_t num, void* h)
  {
    void * memory_block;
    char* ptr = reinterpret_cast<char *>(data) - sizeof(char*);
    size_t taille = sizeof(char*);
    memcpy(&memory_block, ptr, taille);
    for (size_t i = 0; i < num; i++)
      data[i].~T();
    
    free(memory_block);
  }


  template <class T>
  void* MallocObject<T>::reallocate(pointer data, size_t num, void* h)
  {
    if (data == NULL)
      return allocate(num, h);

    void * memory_block; size_t taille = sizeof(char*);
    memcpy(&memory_block,
           reinterpret_cast<char *>(data) - sizeof(char*), taille);

    size_t initial_num = *reinterpret_cast<long*>(memory_block);

    if (initial_num < num)
      {
	memory_block = realloc(memory_block, sizeof(size_t) + sizeof(char*) +
                               (num + 2) * sizeof(T));

	new(static_cast<char *>(memory_block) + sizeof(size_t) + sizeof(T) +
            sizeof(char*) + initial_num * sizeof(T)) T[num - initial_num];
      }
    else if (initial_num > num)
      {
	for (size_t i = num; i < initial_num; i++)
	  data[i].~T();

	memory_block = realloc(memory_block, sizeof(size_t) + sizeof(char*) +
                               (num + 2) * sizeof(T));

      }
    else
      return data;

    memcpy(memory_block, &num, sizeof(size_t));

    pointer data_P =
      reinterpret_cast<pointer>(static_cast<char*>(memory_block) +
                                sizeof(size_t) + sizeof(char*) + sizeof(T));
    memcpy(reinterpret_cast<char *>(data_P) - sizeof(char*),
           &memory_block, taille);

    return data_P;
  }


  template <class T>
  void MallocObject<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }


  template <class T>
  void
  MallocObject<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    for (size_t i = 0; i < num; i++)
      datat[i] = datas[i];
  }


  //////////////
  // NANALLOC //
  //////////////


  template <class T>
  typename NaNAlloc<T>::pointer
  NaNAlloc<T>::allocate(size_t num, void* h)
  {
    pointer data = static_cast<pointer>( malloc(num * sizeof(T)) );
    if (numeric_limits<value_type>::has_signaling_NaN)
      for (size_t i = 0; i < num; i++)
	data[i] = numeric_limits<value_type>::signaling_NaN();
    else if (numeric_limits<value_type>::has_quiet_NaN)
      for (size_t i = 0; i < num; i++)
	data[i] = numeric_limits<value_type>::quiet_NaN();
    else if  (numeric_limits<value_type>::has_infinity)
      for (size_t i = 0; i < num; i++)
	data[i] = numeric_limits<value_type>::infinity();
    return data;
  }

  template <class T>
  void* NaNAlloc<T>::reallocate(pointer data, size_t num, void* h)
  {
    void* datav = realloc(reinterpret_cast<void*>(data), num * sizeof(T));
    pointer datap = reinterpret_cast<pointer>(datav);
    if (numeric_limits<value_type>::has_signaling_NaN)
      for (size_t i = 0; i < num; i++)
	datap[i] = numeric_limits<value_type>::signaling_NaN();
    else if (numeric_limits<value_type>::has_quiet_NaN)
      for (size_t i = 0; i < num; i++)
	datap[i] = numeric_limits<value_type>::quiet_NaN();
    else if  (numeric_limits<value_type>::has_infinity)
      for (size_t i = 0; i < num; i++)
	datap[i] = numeric_limits<value_type>::infinity();
    return datav;
  }

} // namespace Seldon.

#define SELDON_FILE_ALLOCATOR_CXX
#endif
