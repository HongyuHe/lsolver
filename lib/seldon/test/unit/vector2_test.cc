#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<Real_wp> VectReal_wp;

Real_wp threshold;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}


template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > eps) || isnan(abs(x(i)-y(i))))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, class T2, class Alloc, class Alloc2>
bool EqualVector(const Vector<Vector<T>, VectFull, Alloc>& x,
		 const Vector<Vector<T2>, VectFull, Alloc2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if (!EqualVector(x(i), y(i), eps))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, class T2>
bool EqualVector2(const Vector2<T>& x, const Vector2<T2>& y, Real_wp eps = threshold)
{
  if (x.GetLength() != y.GetLength())
    return false;
  
  for (int i = 0; i < x.GetSize(); i++)
    if (!EqualVector(x(i), y(i), eps))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T, class T2>
bool EqualVector3(const Vector3<T>& x, const Vector3<T2>& y, Real_wp eps = threshold)
{
  if (x.GetLength() != y.GetLength())
    return false;
  
  for (int i = 0; i < x.GetSize(); i++)
    if (!EqualVector2(x(i), y(i), eps))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T>
void CheckVector2(const Vector2<T>&)
{
  Vector2<T> x;

  if (!x.IsEmpty())
    {
      cout << "IsEmpty incorrect" << endl;
      abort();
    }
  
  int m = 10;
  Vector<long> length(m);
  for (int i = 0; i < m; i ++)
    length(i) = rand()%(2*m);
  
  x.Reallocate(m);

  if (!x.IsEmpty())
    {
      cout << "IsEmpty incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < m; i ++)
    x.Reallocate(i, length(i));

  if ((x.GetLength() != m) || (x.GetSize() != m))
    {
      cout << "GetLength or GetSize incorrect" << endl;
      abort();
    }

  int sum = 0;
  for (int i = 0; i < m; i ++)
    {
      sum += x.GetLength(i);
      if ((x.GetLength(i) != length(i)) || (x.GetSize(i) != length(i)))
	{
	  cout << "GetLength or GetSize incorrect" << endl;
	  abort();
	}
    }

  if (x.GetMemorySize() < m*sizeof(VectReal_wp) + sizeof(T)*sum)
    {
      cout << "GetMemorySize incorrect" << endl;
      DISP(x.GetMemorySize()); DISP(m*sizeof(VectReal_wp) + sizeof(T)*sum);
      abort();
    }

  if (x.GetNelement() != sum)
    {
      cout << "GetNelement incorrect" << endl;
      abort();
    }

  Vector<long> length2 = x.GetShape();
  if (!EqualVector(length, length2))
    {
      cout << "GetShape incorrect" << endl;
      abort();
    }

  length2.Clear();
  x.GetShape(length2);
  if (!EqualVector(length, length2))
    {
      cout << "GetShape incorrect" << endl;
      abort();
    }

  x.Clear();
  if (!x.IsEmpty() || x.GetSize() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  x.Reallocate(length);
  for (int i = 0; i < m; i ++)
    if (x.GetLength(i) != length(i))
      {
	cout << "Reallocate incorrect" << endl;
	abort();
      }

  Vector<Vector<T> > y(x.GetSize());
  for (int i = 0; i < m; i++)
    {
      y(i).Reallocate(length(i));
      for (int j = 0; j < length(i); j++)
	{
	  GetRandNumber(y(i)(j));
	  x(i, j) = y(i)(j);
	}
    }

  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < length(i); j++)
	if (x(i, j) != y(i)(j))
	  {
	    cout << "Operator () incorrect" << endl;
	    abort();
	  }

      if (!EqualVector(x(i), y(i)))
	{
	  cout << "Operator () incorrect" << endl;
	  abort();
	}

      if (!EqualVector(x.GetVector(i), y(i)))
	{
	  cout << "Operator () incorrect" << endl;
	  abort();
	}
    }

  Vector2<T> z = x.Copy();

  for (int i = 0; i < m; i++)
    {
      if (!EqualVector(x.GetVector(i), z.GetVector(i)))
	{
	  cout << "Copy incorrect" << endl;
	  abort();
	}
    }

  Vector2<T> w;
  w.Copy(x);  
  if (!EqualVector2(x, w))
    {
      cout << "Copy incorrect" << endl;
      abort();
    }

  int beg = rand()%(m-1);
  int end = beg + 1 + rand()%(m-beg-1);
  z.Reallocate(end-beg);
  for (int i = beg; i < end; i++)
    z(i-beg) = x(i);

  x.Select(beg, end);
  if (!EqualVector2(x, z))
    {
      cout << "Select incorrect" << endl;
      abort();
    }

  Vector<T> v(x.GetNelement());
  int nb = 0;
  for (int i = 0; i < x.GetSize(); i++)
    for (int j = 0; j < x.GetSize(i); j++)
      v(nb++) = x(i, j);

  if (!EqualVector(v, x.Flatten()))
    {
      cout << "Flatten incorrect" << endl;
      abort();
    }

  Vector<T> v2;
  x.Flatten(v2);
  if (!EqualVector(v, v2))
    {
      cout << "Flatten incorrect" << endl;
      abort();
    }

  v2.Clear();
  w.Flatten(beg, end, v2);
  if (!EqualVector(v, v2))
    {
      cout << "Flatten incorrect" << endl;
      abort();
    }

  y.Reallocate(x.GetSize());
  for (int i = 0; i < y.GetM(); i++)
    y(i) = x(i);
  
  int p = rand()%x.GetLength();
  T alpha; GetRandNumber(alpha);
  x.PushBack(p, alpha);
  y(p).PushBack(alpha);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }

  Vector<T> vec;
  GenerateRandomVector(vec, rand()%(2*m));

  x.PushBack(vec);
  y.PushBack(vec);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }

  Vector<Vector<T> > yp(m/2+1);

  z.Reallocate(yp.GetM());
  for (int i = 0; i < yp.GetM(); i++)
    {
      GenerateRandomVector(yp(i), rand()%(2*m));
      z(i) = yp(i);
    }

  x.PushBack(z);
  y.PushBack(yp);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }

  y.Reallocate(w.GetSize());
  for (int i = 0; i < y.GetM(); i++)
    y(i) = w(i);

  y(p).Clear();
  w.Clear(p);
  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), w(i)))
      {
	cout << "Clear incorrect" << endl;
	abort();
      }

  w.Clear();
  if (w.GetNelement() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  z.Fill(alpha);
  vec.Reallocate(z.GetNelement());
  vec.Fill(alpha);

  if (!EqualVector(vec, z.Flatten()))
    {
      cout << "Fill incorrect" << endl;
      abort();
    }

  
}


template<class T>
void CheckVector3(const Vector3<T>&)
{
  Vector3<T> x;

  int m = 10;
  Vector<IVect> length(m);
  for (int i = 0; i < m; i ++)
    {
      int n = rand()%(2*m);
      length(i).Reallocate(n);
      for (int j = 0; j < n; j++)
	length(i)(j) = rand()%(2*n);
    }
  
  x.Reallocate(m);  
  for (int i = 0; i < m; i ++)
    {
      x.Reallocate(i, length(i).GetM());
      for (int j = 0; j < length(i).GetM(); j++)
	x.Reallocate(i, j, length(i)(j));
    }
  
  if ((x.GetLength() != m) || (x.GetSize() != m))
    {
      cout << "GetLength or GetSize incorrect" << endl;
      abort();
    }

  int sum = 0;
  for (int i = 0; i < m; i ++)
    {
      if ((x.GetLength(i) != length(i).GetM()) || (x.GetSize(i) != length(i).GetM()))
	{
	  cout << "GetLength or GetSize incorrect" << endl;
	  abort();
	}

      for (int j = 0; j < length(i).GetM(); j++)
	{
	  sum += x.GetLength(i, j);

	  if ((x.GetLength(i, j) != length(i)(j)) || (x.GetSize(i, j) != length(i)(j)))
	    {
	      cout << "GetLength or GetSize incorrect" << endl;
	      abort();
	    }
	}
    }

  if (x.GetMemorySize() < m*sizeof(VectReal_wp) + sizeof(T)*sum)
    {
      cout << "GetMemorySize incorrect" << endl;
      DISP(x.GetMemorySize()); DISP(m*sizeof(VectReal_wp) + sizeof(T)*sum);
      abort();
    }

  if (x.GetNelement() != sum)
    {
      cout << "GetNelement incorrect" << endl;
      abort();
    }
  
  int p = rand()%x.GetLength();
  Vector<int> length2 = x.GetShape(p);
  if (!EqualVector(length(p), length2))
    {
      cout << "GetShape incorrect" << endl;
      abort();
    }

  length2.Clear();
  x.GetShape(p, length2);
  if (!EqualVector(length(p), length2))
    {
      cout << "GetShape incorrect" << endl;
      abort();
    }

  Vector<Vector<Vector<T> > > y(m);
  for (int i = 0; i < m; i++)
    {
      y(i).Reallocate(length(i).GetM());
      for (int j = 0; j < length(i).GetM(); j++)
	{
	  y(i)(j).Reallocate(length(i)(j));
	  for (int k = 0; k < length(i)(j); k++)
	    {
	      GetRandNumber(y(i)(j)(k));
	      x(i, j, k) = y(i)(j)(k);
	    }
	}
    }

  for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < length(i).GetM(); j++)
	for (int k = 0; k < length(i)(j); k++)
	  if (x(i, j, k) != y(i)(j)(k))
	    {
	      cout << "Operator () incorrect" << endl;
	      abort();
	    }

      for (int j = 0; j < length(i).GetM(); j++)
	{
	  if (!EqualVector(x(i, j), y(i)(j)))
	    {
	      cout << "Operator () incorrect" << endl;
	      abort();
	    }
	  
	  if (!EqualVector(x.GetVector(i, j), y(i)(j)))
	    {
	      cout << "Operator () incorrect" << endl;
	      abort();
	    }
	}

      if (!EqualVector(x(i), y(i)))
	{
	  cout << "Operator () incorrect" << endl;
	  abort();
	}
      
      if (!EqualVector(x.GetVector(i), y(i)))
	{
	  cout << "Operator () incorrect" << endl;
	  abort();
	}
    }

  int beg = rand()%(m-1);
  int end = beg + 1 + rand()%(m-beg-1);
  sum = 0;
  for (int i = beg; i < end; i++)
    for (int j = 0; j < length(i).GetM(); j++)
      sum += length(i)(j); 

  if (sum != x.GetNelement(beg, end))
    {
      cout << "GetNelement incorrect" << endl;
      abort();
    }

  Vector3<T> w(x);
  Vector<T> v(x.GetNelement());
  int nb = 0;
  for (int i = 0; i < x.GetSize(); i++)
    for (int j = 0; j < x.GetSize(i); j++)
      for (int k = 0; k < x.GetSize(i, j); k++)
	v(nb++) = x(i, j, k);
  
  Vector<T> v2;
  x.Flatten(v2);
  if (!EqualVector(v, v2))
    {
      cout << "Flatten incorrect" << endl;
      abort();
    }

  v.Reallocate(x.GetNelement(beg, end));
  nb = 0;
  for (int i = beg; i < end; i++)
    for (int j = 0; j < x.GetSize(i); j++)
      for (int k = 0; k < x.GetSize(i, j); k++)
	v(nb++) = x(i, j, k);
  
  v2.Clear();
  w.Flatten(beg, end, v2);
  if (!EqualVector(v, v2))
    {
      cout << "Flatten incorrect" << endl;
      abort();
    }

  y.Reallocate(x.GetSize());
  for (int i = 0; i < y.GetM(); i++)
    y(i) = x(i);
  
  p = rand()%x.GetLength();
  int q = rand()%x.GetLength(p);
  T alpha; GetRandNumber(alpha);
  x.PushBack(p, q, alpha);
  y(p)(q).PushBack(alpha);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
  
  Vector<T> vec;
  GenerateRandomVector(vec, rand()%(2*m));

  x.PushBack(p, vec);
  y(p).PushBack(vec);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }

  Vector<Vector<T> > yp(m/2+1);
  Vector<Vector<T>, VectFull, MallocObject<Vector<T> > > ypb(m/2+1);
  for (int i = 0; i < yp.GetM(); i++)
    {
      GenerateRandomVector(yp(i), rand()%(2*m));
      ypb(i) = yp(i);
    }
  
  x.PushBack(ypb);
  y.PushBack(yp);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }

  Vector<Vector<Vector<T> > > yp2(m/2+1);
  Vector3<T> ypb2(m/2+1);
  for (int i = 0; i < yp2.GetM(); i++)
    {
      int n = rand()%(2*m);
      yp2(i).Reallocate(n);
      for (int j = 0; j < n; j++)
	GenerateRandomVector(yp2(i)(j), rand()%(2*n));

      ypb2(i) = yp2(i);
    }
  
  x.PushBack(ypb2);
  y.PushBack(yp2);

  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), x(i)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }

  y.Reallocate(w.GetSize());
  for (int i = 0; i < y.GetM(); i++)
    y(i) = w(i);

  y(p)(q).Clear();
  w.Clear(p, q);
  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), w(i)))
      {
	cout << "Clear incorrect" << endl;
	abort();
      }
  
  y(p).Clear();
  w.Clear(p);
  for (int i = 0; i < y.GetM(); i++)
    if (!EqualVector(y(i), w(i)))
      {
	cout << "Clear incorrect" << endl;
	abort();
      }
  
  w.Clear();
  if (w.GetNelement() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  x.Fill(alpha);
  vec.Reallocate(x.GetNelement());
  vec.Fill(alpha);

  x.Flatten(v);
  if (!EqualVector(vec, v))
    {
      cout << "Fill incorrect" << endl;
      abort();
    }

  x.Clear();
  if (x.GetSize() != 0)
    {
      cout << "Clear incorrect" << endl;
      abort();
    }  
}

int main(int argc, char** argv)
{
  threshold = 1e-12;

  {
    Vector2<Real_wp> x;
    CheckVector2(x);

    Vector3<Real_wp> x3;
    CheckVector3(x3);
  }

  {
    Vector2<Complex_wp> x;
    CheckVector2(x);
  }

  cout << "All tests passed successfully" << endl;
  
  return 0;
}
