// Copyright (C) 2001-2015 Marc Duruflé
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

#ifndef SELDON_FILE_DEFAULT_ALLOCATOR_HXX

// default allocators for basic types

// if you want to redefine you own default allocators, you should
// include a file with the SELDON_FILE_DEFAULT_ALLOCATOR_HXX
// such that this file is not included

namespace Seldon
{
  template<class Storage>
  class SeldonDefaultAllocator<Storage, bool>
  {
  public:
    typedef MallocAlloc<bool> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, char>
  {
  public:
    typedef MallocAlloc<char> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, unsigned>
  {
  public:
    typedef MallocAlloc<unsigned> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, int>
  {
  public:
    typedef MallocAlloc<int> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, unsigned long>
  {
  public:
    typedef MallocAlloc<unsigned long> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, long>
  {
  public:
    typedef MallocAlloc<long> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, unsigned long long>
  {
  public:
    typedef MallocAlloc<unsigned long long> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, long long>
  {
  public:
    typedef MallocAlloc<long long> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, float>
  {
  public:
    typedef MallocAlloc<float> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, double>
  {
  public:
    typedef MallocAlloc<double> allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, long double>
  {
  public:
    typedef MallocAlloc<long double> allocator;
  };

  template<class Storage>
  class SeldonDefaultAllocator<Storage, complex<float> >
  {
  public:
    typedef MallocAlloc<complex<float> > allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, complex<double> >
  {
  public:
    typedef MallocAlloc<complex<double> > allocator;
  };  

  template<class Storage>
  class SeldonDefaultAllocator<Storage, complex<long double> >
  {
  public:
    typedef MallocAlloc<complex<long double> > allocator;    
  };  
  
}

#define SELDON_FILE_DEFAULT_ALLOCATOR_HXX
#endif
