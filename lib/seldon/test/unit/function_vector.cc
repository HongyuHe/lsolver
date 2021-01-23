#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<Real_wp> VectReal_wp;

namespace std
{
  inline bool isnan(const complex<double>& x)
  {
    if (isnan(real(x)))
      return true;
    
    if (isnan(imag(x)))
      return true;
    
    return false;
  }
}

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

template<class T>
void GenerateRandomVector(Vector<T, VectSparse>& x, int n, int nnz)
{
  x.Clear();
  for (int i = 0; i < nnz; i++)
    {
      int j = rand()%n;
      GetRandNumber(x.Get(j));
    }
}

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > threshold) || isnan(abs(x(i)-y(i))))
      {
        DISP(i); DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T>
bool EqualVector(const Vector<T, VectSparse>& x, const Vector<T, VectSparse>& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if ( (x.Index(i) != y.Index(i)) || (abs(x.Value(i) - y.Value(i)) > threshold) 
         || isnan(abs(x.Value(i) - y.Value(i))) )
      {
        DISP(i); DISP(x.Value(i)); DISP(y.Value(i));
        return false;
      }
  
  return true;
}

int main(int argc, char** argv)
{
  threshold = 2e-14;
  
  //srand(time(NULL));

  {
    // testing real vectors
    int n = 10;
    Vector<Real_wp> x(n), y(n);    
    GenerateRandomVector(x, n);
    Copy(x, y);
    
    if (!EqualVector(x, y))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Real_wp beta = Real_wp(-22)/Real_wp(10);
    Real_wp alpha(beta);
    Mlt(alpha, x);
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(x(i) - beta*y(i)) > threshold) || isnan(abs(x(i)-beta*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    
    Real_wp gamma = -1.2;
    x = y;
    Mlt(gamma, x);
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(x(i) - gamma*y(i)) > threshold) || isnan(abs(x(i) - gamma*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}

    GenerateRandomVector(y, n);
    Vector<Real_wp> z(y);
    Add(alpha, x, y);
    
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(y(i) - (z(i) + alpha*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + alpha*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    y = z;
    Add(gamma, x, y);
    
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(y(i) - (z(i) + gamma*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + gamma*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    z = y;
    Vector<Real_wp> w(x);
    Swap(x, y);
    
    if ( !EqualVector(x, z) || !EqualVector(y, w) )
      {
	cout << "Swap incorrect" << endl;
	abort();
      }
    
    Real_wp scal = DotProd(x, y);
    Real_wp val_ref, zero;
    SetComplexZero(zero);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += x(i)*y(i);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "DotProd incorrect" << endl;
	abort();
      }
    
    scal = Norm1(y);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += ComplexAbs(y(i));
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "Norm1 incorrect" << endl;
	abort();
      }
    
    scal = Norm2(x);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += x(i)*x(i);
    
    val_ref = sqrt(val_ref);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "Norm2 incorrect" << endl;
	abort();
      }    
    
    int imax = GetMaxAbsIndex(x);
    int iref = 0;
    Real_wp xmax = ComplexAbs(x(0));
    for (int i = 1; i < n; i++)
      {
	if (ComplexAbs(x(i)) > xmax)
	  {
	    iref = i;
	    xmax = ComplexAbs(x(i));
	  }
      }
    
    if (iref != imax)
      {
        DISP(x); DISP(iref); DISP(imax);
	cout << "GetMaxAbsIndex incorrect" << endl;
	abort();
      }
    
    Real_wp a(2.1), b(0.7), c, s, a2 = a, b2 = b;
    GenRot(a, b, c, s);
    a = a2; b = b2;
    ApplyRot(a2, b2, c, s);
    if (( (abs(a2 - sqrt(a*a + b*b)) > threshold) || abs(b2) > threshold)
        || isnan(abs(a2 - sqrt(a*a + b*b))) || isnan(abs(b2)))
      {
	cout << "GenRot/ApplyRot incorrect" << endl;
	abort();
      }
      
  }
 
 {
    // testing complex vectors
    int n = 10;
    Vector<Complex_wp> x(n), y(n);    
    GenerateRandomVector(x, n);
    Copy(x, y);
    
    if (!EqualVector(x, y))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Complex_wp alphac(1.9, -0.7);
    Mlt(alphac, x);
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(x(i) - alphac*y(i)) > threshold) || isnan(abs(x(i) - alphac*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    
    Real_wp beta = Real_wp(-22)/Real_wp(10);
    Real_wp alpha(beta);
    x = y;
    Mlt(alpha, x);
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(x(i) - beta*y(i)) > threshold) || isnan(abs(x(i) - beta*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    
    Real_wp gamma = -1.2;
    x = y;
    Mlt(gamma, x);
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(x(i) - Real_wp(gamma)*y(i)) > threshold) || isnan(abs(x(i) - Real_wp(gamma)*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}

    GenerateRandomVector(y, n);
    Vector<Complex_wp> z(y);
    Add(alphac, x, y);
    
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(y(i) - (z(i) + alphac*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + alphac*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    y = z;
    Add(alpha, x, y);
    
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(y(i) - (z(i) + alpha*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + alpha*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    y = z;
    Add(Real_wp(gamma), x, y);
    
    for (int i = 0; i < x.GetM(); i++)
      if ((abs(y(i) - (z(i) + Real_wp(gamma)*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + Real_wp(gamma)*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}

    z = y;
    Vector<Complex_wp> w(x);
    Swap(x, y);
    
    if ( !EqualVector(x, z) || !EqualVector(y, w) )
      {
	cout << "Swap incorrect" << endl;
	abort();
      }
    
    Complex_wp scal = DotProd(x, y);
    Complex_wp val_ref, zero;
    SetComplexZero(zero);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += x(i)*y(i);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
        DISP(val_ref); DISP(scal);
	cout << "DotProd incorrect" << endl;
	abort();
      }
    
    scal = DotProdConj(x, y);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += conjugate(x(i))*y(i);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "DotProdConj incorrect" << endl;
	abort();
      }
    
    Real_wp scal_r = Norm1(y);
    Real_wp val_real(0);
    for (int i = 0; i < n; i++)
      val_real += ComplexAbs(y(i));
    
    if ((abs(val_real - scal_r) > threshold) || isnan(abs(val_real - scal_r)))
      {
	cout << "Norm1 incorrect" << endl;
	abort();
      }
    
    scal_r = Norm2(x);
    val_real = Real_wp(0);
    for (int i = 0; i < n; i++)
      val_real += real(conjugate(x(i))*x(i));
    
    val_real = sqrt(val_real);
    
    if ((abs(val_real - scal_r) > threshold) || isnan(abs(val_real - scal_r)))
      {
	cout << "Norm2 incorrect" << endl;
	abort();
      }    
    
    w = x;
    Conjugate(x);
    for (int i = 0; i < n; i++)
      if ((x(i) != conjugate(w(i))) || isnan(x(i)) || isnan(w(i)))
	{
	  cout << "Conjugate incorrect" << endl;
	  abort();
	}
    
    int imax = GetMaxAbsIndex(x);
    int iref = 0;
    Real_wp xmax = ComplexAbs(x(0));
    for (int i = 1; i < n; i++)
      {
	if (ComplexAbs(x(i)) > xmax)
	  {
	    iref = i;
	    xmax = ComplexAbs(x(i));
	  }
      }
    
    if (iref != imax)
      {
        DISP(x); DISP(iref); DISP(imax);
	cout << "GetMaxAbsIndex incorrect" << endl;
	abort();
      }
    
    Complex_wp a(2.1, 1.3), b(0.7, 1.2), s, a2 = a, b2 = b;
    Complex_wp sqrt_a;
    sqrt_a = sqrt(a);
    Real_wp c;
    GenRot(a, b, c, s);
    a = a2; b = b2;
    ApplyRot(a2, b2, c, s);
    if (( (abs(abs(a2) - sqrt(conjugate(a)*a + conjugate(b)*b)) > threshold) || abs(b2) > threshold)
        || isnan(abs(abs(a2) - sqrt(conjugate(a)*a + conjugate(b)*b))) || isnan(abs(b2)))
      {
        DISP(abs(abs(a2) - sqrt(conjugate(a)*a + conjugate(b)*b)));
        DISP(abs(b2));
	cout << "GenRot/ApplyRot incorrect" << endl;
	abort();
      }
      
  }

  {
    // testing real sparse vectors
    int n = 30, nnz = 10;
    Vector<Real_wp, VectSparse> x, y;    
    GenerateRandomVector(x, n, nnz);
    Copy(x, y);
    
    if (!EqualVector(x, y))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Real_wp beta = Real_wp(-22)/Real_wp(10);
    Real_wp alpha(beta);
    Mlt(alpha, x);
    for (int i = 0; i < n; i++)
      if ((abs(x(i) - beta*y(i)) > threshold) || isnan(abs(x(i) - beta*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    
    Real_wp gamma = -1.2;
    x = y;
    Mlt(gamma, x);
    for (int i = 0; i < n; i++)
      if ((abs(x(i) - Real_wp(gamma)*y(i)) > threshold) || isnan(abs(x(i) - Real_wp(gamma)*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}

    GenerateRandomVector(y, n, nnz-2);
    Vector<Real_wp, VectSparse> z(y);
    Add(alpha, x, y);
    
    for (int i = 0; i < n; i++)
      if ((abs(y(i) - (z(i) + alpha*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + alpha*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    y = z;
    Add(gamma, x, y);
    
    for (int i = 0; i < n; i++)
      if ((abs(y(i) - (z(i) + Real_wp(gamma)*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + Real_wp(gamma)*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}

    z = y;
    Vector<Real_wp, VectSparse> w(x);
    Swap(x, y);
    
    if ( !EqualVector(x, z) || !EqualVector(y, w) )
      {
	cout << "Swap incorrect" << endl;
	abort();
      }
    
    Real_wp scal = DotProd(x, y);
    Real_wp val_ref, zero;
    SetComplexZero(zero);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += x(i)*y(i);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "DotProd incorrect" << endl;
	abort();
      }
    
    scal = Norm1(y);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += abs(y(i));
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "Norm1 incorrect" << endl;
	abort();
      }
    
    scal = Norm2(x);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += x(i)*x(i);
    
    val_ref = sqrt(val_ref);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "Norm2 incorrect" << endl;
	abort();
      }    
    
    int imax = GetMaxAbsIndex(x);
    int iref = 0;
    Real_wp xmax = abs(x(0));
    for (int i = 1; i < n; i++)
      {
	if (abs(x(i)) > xmax)
	  {
	    iref = i;
	    xmax = abs(x(i));
	  }
      }
    
    if (iref != imax)
      {
	cout << "GetMaxAbsIndex incorrect" << endl;
	abort();
      }
  }
  
  {
    // testing complex sparse vectors
    int n = 30, nnz = 10;
    Vector<Complex_wp, VectSparse> x, y;    
    GenerateRandomVector(x, n, nnz);
    Copy(x, y);
    
    if (!EqualVector(x, y))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Complex_wp alphac(1.9, -0.7);
    Mlt(alphac, x);
    for (int i = 0; i < n; i++)
      if ((abs(x(i) - alphac*y(i)) > threshold) || isnan(abs(x(i) - alphac*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    
    Real_wp beta = Real_wp(-22)/Real_wp(10);
    Real_wp alpha(beta);
    x = y;
    Mlt(alpha, x);
    for (int i = 0; i < n; i++)
      if ((abs(x(i) - beta*y(i)) > threshold) || isnan(abs(x(i) - beta*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
    
    Real_wp gamma = -1.2;
    x = y;
    Mlt(gamma, x);
    for (int i = 0; i < n; i++)
      if ((abs(x(i) - Real_wp(gamma)*y(i)) > threshold) || isnan(abs(x(i) - Real_wp(gamma)*y(i))))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}

    GenerateRandomVector(y, n, nnz-2);
    Vector<Complex_wp, VectSparse> z(y);
    Add(alphac, x, y);
    
    for (int i = 0; i < n; i++)
      if ((abs(y(i) - (z(i) + alphac*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + alphac*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    y = z;
    Add(alpha, x, y);
    
    for (int i = 0; i < n; i++)
      if ((abs(y(i) - (z(i) + alpha*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + alpha*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
    
    y = z;
    Add(gamma, x, y);
    
    for (int i = 0; i < n; i++)
      if ((abs(y(i) - (z(i) + Real_wp(gamma)*x(i))) > threshold) || isnan(abs(y(i) - (z(i) + Real_wp(gamma)*x(i)))))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}

    z = y;
    Vector<Complex_wp, VectSparse> w(x);
    Swap(x, y);
    
    if ( !EqualVector(x, z) || !EqualVector(y, w) )
      {
	cout << "Swap incorrect" << endl;
	abort();
      }
    
    Complex_wp scal = DotProd(x, y);
    Complex_wp val_ref, zero;
    SetComplexZero(zero);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += x(i)*y(i);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal)))
      {
	cout << "DotProd incorrect" << endl;
	abort();
      }
    
    scal = DotProdConj(x, y);
    val_ref = zero;
    for (int i = 0; i < n; i++)
      val_ref += conjugate(x(i))*y(i);
    
    if ((abs(val_ref - scal) > threshold) || isnan(abs(val_ref - scal) ))
      {
	cout << "DotProdConj incorrect" << endl;
	abort();
      }
    
    Real_wp scal_r = Norm1(y);
    Real_wp val_real(0);
    for (int i = 0; i < n; i++)
      val_real += ComplexAbs(y(i));
    
    if ((abs(val_real - scal_r) > threshold) || isnan(abs(val_real - scal_r)))
      {
	cout << "Norm1 incorrect" << endl;
	abort();
      }
    
    scal_r = Norm2(x);
    val_real = Real_wp(0);
    for (int i = 0; i < n; i++)
      val_real += real(conjugate(x(i))*x(i));
    
    val_real = sqrt(val_real);
    
    if ((abs(val_real - scal_r) > threshold) || isnan(abs(val_real - scal_r)))
      {
	cout << "Norm2 incorrect" << endl;
	abort();
      }    

    w = x;
    Conjugate(x);
    for (int i = 0; i < n; i++)
      if ((x(i) != conjugate(w(i))) || isnan(x(i)) || isnan(w(i)))
	{
	  cout << "Conjugate incorrect" << endl;
	  abort();
	}
    
    int imax = GetMaxAbsIndex(x);
    int iref = 0;
    Real_wp xmax = ComplexAbs(x(0));
    for (int i = 1; i < n; i++)
      {
	if (ComplexAbs(x(i)) > xmax)
	  {
	    iref = i;
	    xmax = abs(x(i));
	  }
      }
    
    if (iref != imax)
      {
	cout << "GetMaxAbsIndex incorrect" << endl;
	abort();
      }
    
   }
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
