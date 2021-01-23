// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_ARRAYS_HXX

namespace Seldon
{

  template<class T, class Storage, class Allocator>
  long PartitionQuickSort(long m, long n,
			 Vector<T, Storage, Allocator>& t);

  template<class T, class Storage, class Allocator>
  void QuickSort(long m, long n,
		 Vector<T, Storage, Allocator>& t);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  long PartitionQuickSort(long m, long n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void QuickSort(long m, long n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  long PartitionQuickSort(long m, long n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2,
			 Vector<T3, Storage3, Allocator3>& t3);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void QuickSort(long m, long n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2,
		 Vector<T3, Storage3, Allocator3>& t3);
  
  template<class T, class Storage, class Allocator>
  void MergeSort(long m, long n, Vector<T, Storage, Allocator>& tab1);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void MergeSort(long m, long n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void MergeSort(long m, long n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2,
		 Vector<T3, Storage3, Allocator3>& tab3);
  
  template<class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2 >
  void Assemble(long& n, Vector<int, Storage1, Allocator1>& Node,
		Vector<T2, Storage2, Allocator2>& Vect);
  
  template<class T, class Storage1, class Allocator1>
  void Assemble(long& n, Vector<T, Storage1, Allocator1>& Node);
  
  template<class T, class Storage1, class Allocator1>
  void Assemble(Vector<T, Storage1, Allocator1>& Node);

  template<class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2 >
  inline void Assemble(int& n, Vector<int, Storage1, Allocator1>& Node,
                       Vector<T2, Storage2, Allocator2>& Vect)
  { long n2(n); Assemble(n2, Node, Vect); n = n2; }
  
  template<class T, class Storage1, class Allocator1>
  inline void Assemble(int& n, Vector<T, Storage1, Allocator1>& Node)
  { long n2(n); Assemble(n2, Node); n = n2; }

  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(long& n, Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2);

  template<class T, class Storage1, class Allocator1>
  void RemoveDuplicate(long& n, Vector<T, Storage1, Allocator1>& Node);

  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2);

  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  inline void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2)
  { long n2(n); RemoveDuplicate(n2, Node, Node2); n = n2; }

  template<class T, class Storage1, class Allocator1>
  inline void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node)
  { long n2(n); RemoveDuplicate(n2, Node); n = n2; }

  template<class T, class Storage1, class Allocator1>
  void RemoveDuplicate(Vector<T, Storage1, Allocator1>& Node);

  template<class T, class Storage, class Allocator>
  void Sort(long m, long n, Vector<T, Storage, Allocator>& V);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(long m, long n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(long m, long n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3);

  template<class T, class Storage, class Allocator>
  void Sort(long n, Vector<T, Storage, Allocator>& V);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(long n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(long n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3);
  
  template<class T, class Storage, class Allocator>
  void Sort(Vector<T, Storage, Allocator>& V);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3);

  template<class T, class Storage, class Allocator>
  bool HasElement(const Vector<T, Storage, Allocator>& X, const T& a);
  
} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_ARRAYS_HXX
#endif
