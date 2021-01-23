// Copyright (C) 2001-2009 Vivien Mallet, Marc Fragu
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


#ifndef SELDON_FILE_ALLOCATOR_HXX

namespace Seldon
{


  /////////////////
  // MALLOCALLOC //
  /////////////////


  template <class T>
  class MallocAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    
    //! true if the allocator keeps previous elements when
    //! calling reallocate
    static const bool KeepDataReallocate = true;

  public:

    static pointer allocate(size_t num, void* h = 0);
    static void deallocate(pointer data, size_t num, void* h = 0);
    static void* reallocate(pointer data, size_t num, void* h = 0);
    static void memoryset(pointer data, char c, size_t num);
    static void memorycpy(pointer datat, pointer datas, size_t num);
    
  };


  /////////////////
  // CALLOCALLOC //
  /////////////////


  template <class T>
  class CallocAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

    //! true if the allocator keeps previous elements when
    //! calling reallocate
    static const bool KeepDataReallocate = true;

  public:

    static pointer allocate(size_t num, void* h = 0);
    static void deallocate(pointer data, size_t num, void* h = 0);
    static void* reallocate(pointer data, size_t num, void* h = 0);
    static void memoryset(pointer data, char c, size_t num);
    static void memorycpy(pointer datat, pointer datas, size_t num);
  };


  //////////////
  // NEWALLOC //
  //////////////


  template <class T>
  class NewAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

    //! true if the allocator keeps previous elements when
    //! calling reallocate
    static const bool KeepDataReallocate = false;

  public:

    static pointer allocate(size_t num, void* h = 0);
    static void deallocate(pointer data, size_t num, void* h = 0);
    static void* reallocate(pointer data, size_t num, void* h = 0);
    static void memoryset(pointer data, char c, size_t num);
    static void memorycpy(pointer datat, pointer datas, size_t num);
  };
  

  //////////////////
  // MallocObject //
  //////////////////


  template <class T>
  class MallocObject
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

    //! true if the allocator keeps previous elements when
    //! calling reallocate
    static const bool KeepDataReallocate = false;
    
  public:

    static pointer allocate(size_t num, void* h = 0);
    static void deallocate(pointer data, size_t num, void* h = 0);
    static void* reallocate(pointer data, size_t num, void* h = 0);
    static void memoryset(pointer data, char c, size_t num);
    static void memorycpy(pointer datat, pointer datas, size_t num);
  };
  
  
  //////////////
  // NANALLOC //
  //////////////


  template <class T>
  class NaNAlloc
  {
  public:
    typedef T value_type;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef T& reference;
    typedef const T& const_reference;

    //! true if the allocator keeps previous elements when
    //! calling reallocate
    static const bool KeepDataReallocate = true;

  public:

    static pointer allocate(size_t num, void* h = 0);
    static void deallocate(pointer data, size_t num, void* h = 0);
    static void* reallocate(pointer data, size_t num, void* h = 0);
    static void memoryset(pointer data, char c, size_t num);
    static void memorycpy(pointer datat, pointer datas, size_t num);
  };

  
  //! Selection of default allocator depending on storage and value type
  template<class Storage, class T>
  class SeldonDefaultAllocator
  {
  public :
    // default choice : SELDON_DEFAULT_ALLOCATOR
    typedef SELDON_DEFAULT_ALLOCATOR<T> allocator;
    
  };

} // namespace Seldon.

#define SELDON_FILE_ALLOCATOR_HXX
#endif
