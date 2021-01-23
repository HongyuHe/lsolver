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


template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > eps) || isnan(abs(x(i)-y(i))))
      {
        DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T>
void CheckOperator(Vector<T>& u)
{
  int n = 10;
  u.Reallocate(n);
  Vector<T> v(n), w(n), res(n);
  GenerateRandomVector(u, n);
  GenerateRandomVector(v, n);
  GenerateRandomVector(w, n);
  
  T alpha, beta;
  GetRandNumber(alpha); GetRandNumber(beta);

  for (int i = 0; i < u.GetM(); i++)
    res(i) = alpha*u(i);
  
  u *= alpha;

  if (!EqualVector(u, res))
    {
      cout << "Operator *= incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = v(i) + w(i);

  u = v + w;

  if (!EqualVector(u, res))
    {
      cout << "Operator + incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = v(i) - w(i);

  u = v - w;

  if (!EqualVector(u, res))
    {
      cout << "Operator - incorrect" << endl;
      abort();
    }


  for (int i = 0; i < u.GetM(); i++)
    res(i) = alpha*v(i) + w(i)*beta;

  u = alpha*v + w*beta;

  if (!EqualVector(u, res))
    {
      cout << "Operator = incorrect" << endl;
      abort();
    }

  GenerateRandomVector(v, n);
  GenerateRandomVector(w, n);

  for (int i = 0; i < u.GetM(); i++)
    res(i) += v(i) - alpha*w(i);

  u += v - alpha*w;

  if (!EqualVector(u, res))
    {
      cout << "Operator += incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) -= v(i) - beta*(w(i) + v(i));

  u -= v - beta*(w+v);

  if (!EqualVector(u, res))
    {
      cout << "Operator -= incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = -v(i);

  u = -v;

  if (!EqualVector(u, res))
    {
      cout << "Operator - incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = alpha  + v(i);

  u = alpha + v;

  if (!EqualVector(u, res))
    {
      cout << "Operator + incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = alpha  + v(i);

  u = v + alpha;

  if (!EqualVector(u, res))
    {
      cout << "Operator + incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = alpha  - v(i);

  u = alpha - v;

  if (!EqualVector(u, res))
    {
      cout << "Operator - incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = -alpha  + v(i);

  u = v - alpha;

  if (!EqualVector(u, res))
    {
      cout << "Operator - incorrect" << endl;
      abort();
    }

  typename ClassComplexType<T>::Treal alpha_real;
  GetRandNumber(alpha_real);
  for (int i = 0; i < u.GetM(); i++)
    res(i) += v(i) - alpha_real*w(i);

  u += v - alpha_real*w;

  if (!EqualVector(u, res))
    {
      cout << "Operator * incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) += v(i) - alpha_real*w(i);

  u += v - w*alpha_real;

  if (!EqualVector(u, res))
    {
      cout << "Operator * incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = v(i)*w(i);

  u = v*w;

  if (!EqualVector(u, res))
    {
      cout << "Operator * incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = v(i)/w(i);

  u = v/w;

  if (!EqualVector(u, res))
    {
      cout << "Operator / incorrect" << endl;
      abort();
    }


  for (int i = 0; i < u.GetM(); i++)
    res(i) += v(i) - alpha_real / w(i);

  u += v - alpha_real / w;

  if (!EqualVector(u, res))
    {
      cout << "Operator / incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) += v(i) - w(i) / alpha_real;

  u += v - w/alpha_real;

  if (!EqualVector(u, res))
    {
      cout << "Operator / incorrect" << endl;
      abort();
    }


  for (int i = 0; i < u.GetM(); i++)
    res(i) = v(i) - alpha / w(i);

  u = v - alpha / w;

  if (!EqualVector(u, res))
    {
      cout << "Operator / incorrect" << endl;
      abort();
    }

  for (int i = 0; i < u.GetM(); i++)
    res(i) = v(i) - w(i) / alpha;

  u = v - w/alpha;

  if (!EqualVector(u, res))
    {
      cout << "Operator / incorrect" << endl;
      abort();
    }

}

int main(int argc, char** argv)
{
  cout.precision(15);

  threshold = 1e-12;

  
  {
    // testing vector of strings
    Vector<string> v;
    if (v.GetM() != 0)
      {
	cout << "non-empty vector" << endl;
	abort();
      }
    
    v.Reallocate(3);
    if ( (v.GetM() != 3) || (v.GetLength() != 3) || (v.GetSize() != 3) )
      {
	cout << "GetM incorrect" << endl;
	abort();
      }
    
    v(0) = string("test.mesh");
    v(1) = string("NONE");
    v(2) = string("TEST_String");
    
    v.Resize(5);
    if ( (v(0) != "test.mesh") || (v(1) != "NONE") || (v(2) != "TEST_String"))
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    v(3) = "socket";
    v(4) = "noname.txt";
    
    v.PushBack(string("toto"));
    if ( (v(0) != "test.mesh") || (v(1) != "NONE") || (v(2) != "TEST_String")
	 || (v(3) != "socket") || (v(4) != "noname.txt") || (v(5) != "toto"))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    v.Resize(4);
    if ( (v(0) != "test.mesh") || (v(1) != "NONE") || (v(2) != "TEST_String")
	 || (v(3) != "socket") )
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    Vector<string> w(2);
    w(0) = "tree"; w(1) = "apple";
    if (w.GetM() != 2)
      {
	cout << "Constructor incorrect" << endl;
	abort();
      }
    
    v.PushBack(w);
    if ( (v(0) != "test.mesh") || (v(1) != "NONE") || (v(2) != "TEST_String")
	 || (v(3) != "socket") || (v(4) != "tree") || (v(5) != "apple") )
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    w.Clear();
    if (w.GetM() != 0)
      {
	cout << "Clear incorrect" << endl;
	abort();
      }
    
    w.SetData(v.GetM(), v.GetData());
    v.Nullify();
    if ( (w(0) != "test.mesh") || (w(1) != "NONE") || (w(2) != "TEST_String")
	 || (w(3) != "socket") || (w(4) != "tree") || (w(5) != "apple") )
      {
	cout << "SetData incorrect" << endl;
	abort();
      }
    
    v = w;
    if ( (v(0) != "test.mesh") || (v(1) != "NONE") || (v(2) != "TEST_String")
	 || (v(3) != "socket") || (v(4) != "tree") || (v(5) != "apple") )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    w.Reallocate(3);
    w.Fill(string("same_value"));
    if ( (w(0) != "same_value") || (w(1) != "same_value") || (w(2) != "same_value") )
      {
	cout << "Fill incorrect" << endl;
	abort();
      }
    
    w.Clear();
    v.WriteText("toto.dat");
    w.ReadText("toto.dat");
    if ( (w(0) != "test.mesh") || (w(1) != "NONE") || (w(2) != "TEST_String")
	 || (w(3) != "socket") || (w(4) != "tree") || (w(5) != "apple") )
      {
	cout << "SetData incorrect" << endl;
	abort();
      }        
  }
  

  {
    // testing vector of Real_wps
    Vector<Real_wp> v;
    if (v.GetM() != 0)
      {
	cout << "non-empty vector" << endl;
	abort();
      }
    
    v.Reallocate(3);
    if ( (v.GetM() != 3) || (v.GetLength() != 3) || (v.GetSize() != 3) )
      {
	cout << "GetM incorrect" << endl;
	abort();
      }
    
    v(0) = 0.25;
    v(1) = 2.5;
    v(2) = -0.125;
    
    v.Resize(5);
    if ( (v(0) != Real_wp(0.25)) || (v(1) != Real_wp(2.5)) || (v(2) != Real_wp(-0.125))
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) )
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    v(3) = 1.25;
    v(4) = 0.75;
    
    v.PushBack(-3.0);
    if ( (v(0) != Real_wp(0.25)) || (v(1) != Real_wp(2.5)) || (v(2) != Real_wp(-0.125))
	 || (v(3) != Real_wp(1.25)) || (v(4) != Real_wp(0.75)) || (v(5) != Real_wp(-3.0))
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) || isnan(v(4)) || isnan(v(5)) )
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    v.Resize(4);
    if ( (v(0) != Real_wp(0.25)) || (v(1) != Real_wp(2.5)) || (v(2) != Real_wp(-0.125))
	 || (v(3) != Real_wp(1.25)) || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) )
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    Vector<Real_wp> w(2);
    w(0) = 2.0625; w(1) = -4.5;
    if (w.GetM() != 2)
      {
	cout << "Constructor incorrect" << endl;
	abort();
      }
    
    v.PushBack(w);
    if ( (v(0) != Real_wp(0.25)) || (v(1) != Real_wp(2.5)) || (v(2) != Real_wp(-0.125))
	 || (v(3) != Real_wp(1.25)) || (v(4) != Real_wp(2.0625)) || (v(5) != Real_wp(-4.5)) 
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) || isnan(v(4)) || isnan(v(5)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    w.Clear();
    if (w.GetM() != 0)
      {
	cout << "Clear incorrect" << endl;
	abort();
      }
    
    w.SetData(v.GetM(), v.GetData());
    v.Nullify();
    if ( (w(0) != Real_wp(0.25)) || (w(1) != Real_wp(2.5)) || (w(2) != Real_wp(-0.125))
	 || (w(3) != Real_wp(1.25)) || (w(4) != Real_wp(2.0625)) || (w(5) != Real_wp(-4.5))
         || isnan(w(0)) || isnan(w(1)) || isnan(w(2)) || isnan(w(3)) || isnan(w(4)) || isnan(w(5)))
      {
	cout << "SetData incorrect" << endl;
	abort();
      }
    
    v = w;
    if ( (v(0) != Real_wp(0.25)) || (v(1) != Real_wp(2.5)) || (v(2) != Real_wp(-0.125))
	 || (v(3) != Real_wp(1.25)) || (v(4) != Real_wp(2.0625)) || (v(5) != Real_wp(-4.5))
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) || isnan(v(4)) || isnan(v(5)))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    w.Reallocate(3);
    w.Fill(-5.5);
    if ( (w(0) != Real_wp(-5.5)) || (w(1) != Real_wp(-5.5)) || (w(2) != Real_wp(-5.5))
         || isnan(w(0)) || isnan(w(1)) || isnan(w(2)))
      {
	cout << "Fill incorrect" << endl;
	abort();
      }
    
    w.Clear();
    v.WriteText("toto.dat");
    w.ReadText("toto.dat");
    if ( (w(0) != Real_wp(0.25)) || (w(1) != Real_wp(2.5)) || (w(2) != Real_wp(-0.125))
	 || (w(3) != Real_wp(1.25)) || (w(4) != Real_wp(2.0625)) || (w(5) != Real_wp(-4.5)) 
         || isnan(w(0)) || isnan(w(1)) || isnan(w(2)) || isnan(w(3)) || isnan(w(4)) || isnan(w(5)))
      {
	cout << "ReadText/WriteText incorrect" << endl;
	abort();
      }        
    
    w.Clear();
    v.Write("totob.dat");
    w.Read("totob.dat");
    if ( (w(0) != 0.25) || (w(1) != 2.5) || (w(2) != -0.125)
	 || (w(3) != 1.25) || (w(4) != 2.0625) || (w(5) != -4.5) )
      {
	cout << "Read/Write incorrect" << endl;
	abort();
      }        
    
    if ((w.GetNormInf() != Real_wp(4.5)) || isnan(w.GetNormInf()))
      {
	cout << "GetNormInf incorrect" << endl;
	abort();
      }

    if (w.GetNormInfIndex() != 5)
      {
	cout << "GetNormInfIndex incorrect" << endl;
	abort();
      }
    
    w.Fill();
    if ((w.GetNormInf() != Real_wp(5.0)) || isnan(w.GetNormInf()))
      {
	cout << "Fill incorrect" << endl;
	abort();
      }
    
    CheckOperator(v);
  }
  
  {
    // testing vector of complex Real_wps
    Vector<complex<Real_wp> > v;
    if (v.GetM() != 0)
      {
	cout << "non-empty vector" << endl;
	abort();
      }
    
    v.Reallocate(3);
    if ( (v.GetM() != 3) || (v.GetLength() != 3) || (v.GetSize() != 3) )
      {
	cout << "GetM incorrect" << endl;
	abort();
      }
    
    complex<Real_wp> a(0.25, 0.5), b(2.5, -3.5), c(-0.125, 1.0);
    v(0) = a;
    v(1) = b;
    v(2) = c;
    
    v.Resize(5);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c))
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    complex<Real_wp> x(1.25, 2.0), y(0.75, 0.0), z(-3.0, 0.125);
    v(3) = x;
    v(4) = y;
    
    v.PushBack(z);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || (v(4) != y) || (v(5) != z)
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) || isnan(v(4)) || isnan(v(5)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    v.Resize(4);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)))
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    y = complex<Real_wp>(2.0625, -1.0);
    z = complex<Real_wp>(-4.5, 2.5);
    Vector<complex<Real_wp> > w(2);
    w(0) = y; w(1) = z;
    if (w.GetM() != 2)
      {
	cout << "Constructor incorrect" << endl;
	abort();
      }
    
    v.PushBack(w);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || (v(4) != y) || (v(5) != z) 
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) || isnan(v(4)) || isnan(v(5)))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    w.Clear();
    if (w.GetM() != 0)
      {
	cout << "Clear incorrect" << endl;
	abort();
      }
    
    w.SetData(v.GetM(), v.GetData());
    v.Nullify();
    if ( (w(0) != a) || (w(1) != b) || (w(2) != c)
	 || (w(3) != x) || (w(4) != y) || (w(5) != z) 
         || isnan(w(0)) || isnan(w(1)) || isnan(w(2)) || isnan(w(3)) || isnan(w(4)) || isnan(w(5)))
      {
	cout << "SetData incorrect" << endl;
	abort();
      }
    
    v = w;
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || (v(4) != y) || (v(5) != z) 
         || isnan(v(0)) || isnan(v(1)) || isnan(v(2)) || isnan(v(3)) || isnan(v(4)) || isnan(v(5)))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    w.Reallocate(3);
    complex<Real_wp> d(-5.5, 1.5);
    w.Fill(d);
    if ( (w(0) != d) || (w(1) != d) || (w(2) != d) || isnan(w(0)) || isnan(w(1)) || isnan(w(2)) )
      {
	cout << "Fill incorrect" << endl;
	abort();
      }
    
    w.Clear();
    v.WriteText("toto.dat");
    w.ReadText("toto.dat");
    if ( (w(0) != a) || (w(1) != b) || (w(2) != c)
	 || (w(3) != x) || (w(4) != y) || (w(5) != z) 
         || isnan(w(0)) || isnan(w(1)) || isnan(w(2)) || isnan(w(3)) || isnan(w(4)) || isnan(w(5)))
      {
	cout << "ReadText/WriteText incorrect" << endl;
	abort();
      }        
    
    w.Clear();
    v.Write("totob.dat");
    w.Read("totob.dat");
    if ( (w(0) != a) || (w(1) != b) || (w(2) != c)
	 || (w(3) != x) || (w(4) != y) || (w(5) != z) )
      {
	cout << "Read/Write incorrect" << endl;
	abort();
      }        
    
    if ((w.GetNormInf() != abs(z)) || isnan(w.GetNormInf()))
      {
	cout << "GetNormInf incorrect" << endl;
	abort();
      }

    if (w.GetNormInfIndex() != 5)
      {
	cout << "GetNormInfIndex incorrect" << endl;
	abort();
      }
    
    w.Fill();
    if ((w.GetNormInf() != Real_wp(5.0)) || isnan(w.GetNormInf()))
      {
	cout << "Fill incorrect" << endl;
	abort();
      }    

    CheckOperator(v);
  }
   
  
  {
    // testing vector of integers
    Vector<int> v;
    if (v.GetM() != 0)
      {
	cout << "non-empty vector" << endl;
	abort();
      }
    
    v.Reallocate(3);
    if ( (v.GetM() != 3) || (v.GetLength() != 3) || (v.GetSize() != 3) )
      {
	cout << "GetM incorrect" << endl;
	abort();
      }
    
    int a(2), b(3), c(-1);
    v(0) = a;
    v(1) = b;
    v(2) = c;
    
    v.Resize(5);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c))
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    int x(1), y(7), z(-3);
    v(3) = x;
    v(4) = y;
    
    v.PushBack(z);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || (v(4) != y) || (v(5) != z))
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    v.Resize(4);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) )
      {
	cout << "Resize incorrect" << endl;
	abort();
      }
    
    y = 6;
    z = 4;
    Vector<int> w(2);
    w(0) = y; w(1) = z;
    if (w.GetM() != 2)
      {
	cout << "Constructor incorrect" << endl;
	abort();
      }
    
    v.PushBack(w);
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || (v(4) != y) || (v(5) != z) )
      {
	cout << "PushBack incorrect" << endl;
	abort();
      }
    
    w.Clear();
    if (w.GetM() != 0)
      {
	cout << "Clear incorrect" << endl;
	abort();
      }
    
    w.SetData(v.GetM(), v.GetData());
    v.Nullify();
    if ( (w(0) != a) || (w(1) != b) || (w(2) != c)
	 || (w(3) != x) || (w(4) != y) || (w(5) != z) )
      {
	cout << "SetData incorrect" << endl;
	abort();
      }
    
    v = w;
    if ( (v(0) != a) || (v(1) != b) || (v(2) != c)
	 || (v(3) != x) || (v(4) != y) || (v(5) != z) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    w.Reallocate(3);
    int d(13);
    w.Fill(d);
    if ( (w(0) != d) || (w(1) != d) || (w(2) != d) )
      {
	cout << "Fill incorrect" << endl;
	abort();
      }
    
    w.Clear();
    v.WriteText("toto.dat");
    w.ReadText("toto.dat");
    if ( (w(0) != a) || (w(1) != b) || (w(2) != c)
	 || (w(3) != x) || (w(4) != y) || (w(5) != z) )
      {
	cout << "ReadText/WriteText incorrect" << endl;
	abort();
      }        
    
    w.Clear();
    v.Write("totob.dat");
    w.Read("totob.dat");
    if ( (w(0) != a) || (w(1) != b) || (w(2) != c)
	 || (w(3) != x) || (w(4) != y) || (w(5) != z) )
      {
	cout << "Read/Write incorrect" << endl;
	abort();
      }        
    
    if (w.GetNormInf() != abs(y))
      {
	cout << "GetNormInf incorrect" << endl;
	abort();
      }

    if (w.GetNormInfIndex() != 4)
      {
	cout << "GetNormInfIndex incorrect" << endl;
	abort();
      }
    
    w.Fill();
    if (w.GetNormInf() != 5)
      {
	cout << "Fill incorrect" << endl;
	abort();
      }    
  }
  
  {
    // testing sparse vector of Real_wps
    Vector<Real_wp, VectSparse> v;
    
    v.Reallocate(4);
    if (v.GetM() != 4)
      {
	cout << "GetM incorrect" << endl;
	abort();
      }
    
    v.Index(0) = 3;
    v.Index(1) = 1;
    v.Index(2) = 5;
    v.Index(3) = 3;
    v.Value(0) = to_num<Real_wp>("1.2");
    v.Value(1) = to_num<Real_wp>("0.8");
    v.Value(2) = to_num<Real_wp>("2.7");
    v.Value(3) = to_num<Real_wp>("-0.3");
    
    v.Assemble();
    if ( (v.GetM() != 3) || (v.Index(0) != 1) || (v.Index(1) != 3)
	 || (v.Index(2) != 5) || (abs(v.Value(0) - to_num<Real_wp>("0.8")) > threshold)
	 || (abs(v.Value(1) - to_num<Real_wp>("0.9")) > threshold) || (abs(v.Value(2) - to_num<Real_wp>("2.7")) > threshold) 
         || isnan(v.Value(0)) || isnan(v.Value(1)) || isnan(v.Value(2)))
      {
        DISP(v.Value(1) - to_num<Real_wp>("0.9"));
        DISP(v.Value(0) - to_num<Real_wp>("0.8")); DISP(v.Value(2) - to_num<Real_wp>("2.7"));
	cout << "Assemble incorrect" << endl;
	abort();
      }
    
    v.Resize(5);
    v.Index(3) = 8;
    v.Value(3) = to_num<Real_wp>("3.1");
    v.Index(4) = 2;
    v.Value(4) = to_num<Real_wp>("-2.6");
    v.Assemble();
    
    if ( (v.GetM() != 5) || (v.Index(0) != 1) || (v.Index(1) != 2) || (v.Index(2) != 3)
	 || (v.Index(3) != 5) || (v.Index(4) != 8)
	 || (abs(v.Value(0) - to_num<Real_wp>("0.8")) > threshold)
	 || (abs(v.Value(1) + to_num<Real_wp>("2.6")) > threshold) || (abs(v.Value(2) - to_num<Real_wp>("0.9")) > threshold)
	 || (abs(v.Value(3) - to_num<Real_wp>("2.7")) > threshold) || (abs(v.Value(4) - to_num<Real_wp>("3.1")) > threshold) 
         || isnan(v.Value(0)) || isnan(v.Value(1)) || isnan(v.Value(2))
         || isnan(v.Value(3)) || isnan(v.Value(4)))
      {
	cout << "Resize/ Assemble incorrect" << endl;
	abort();
      }
    
    Vector<Real_wp, VectSparse> w;
    w.SetData(v);
    v.Nullify();
    
    if ( (w.GetM() != 5) || (w.Index(0) != 1) || (w.Index(1) != 2) || (w.Index(2) != 3)
	 || (w.Index(3) != 5) || (w.Index(4) != 8)
	 || (abs(w.Value(0) - to_num<Real_wp>("0.8")) > threshold)
	 || (abs(w.Value(1) + to_num<Real_wp>("2.6")) > threshold) || (abs(w.Value(2) - to_num<Real_wp>("0.9")) > threshold)
	 || (abs(w.Value(3) - to_num<Real_wp>("2.7")) > threshold) || (abs(w.Value(4) - to_num<Real_wp>("3.1")) > threshold) 
         || isnan(w.Value(0)) || isnan(w.Value(1)) || isnan(w.Value(2))
         || isnan(w.Value(3)) || isnan(w.Value(4)))
      {
	cout << "SetData incorrect" << endl;
	abort();
      }
    
    Vector<Real_wp, VectSparse> u(w);
    
    if ( (u.GetM() != 5) || (u.Index(0) != 1) || (u.Index(1) != 2) || (u.Index(2) != 3)
	 || (u.Index(3) != 5) || (u.Index(4) != 8)
	 || (abs(u.Value(0) - to_num<Real_wp>("0.8")) > threshold)
	 || (abs(u.Value(1) + to_num<Real_wp>("2.6")) > threshold) || (abs(u.Value(2) - to_num<Real_wp>("0.9")) > threshold)
	 || (abs(u.Value(3) - to_num<Real_wp>("2.7")) > threshold) || (abs(u.Value(4) - to_num<Real_wp>("3.1")) > threshold) 
         || isnan(u.Value(0)) || isnan(u.Value(1)) || isnan(u.Value(2))
         || isnan(u.Value(3)) || isnan(u.Value(4)))
      {
	cout << "Constructor incorrect" << endl;
	abort();
      }
    
    u.Get(2) += to_num<Real_wp>("1.4");
    u.Get(7) = to_num<Real_wp>("-2.9");
    
    if ( (u.GetM() != 6) || (u.Index(0) != 1) || (u.Index(1) != 2) || (u.Index(2) != 3)
	 || (u.Index(3) != 5) || (u.Index(4) != 7) || (u.Index(5) != 8)
	 || (abs(u.Value(0) - to_num<Real_wp>("0.8")) > threshold) || (abs(u.Value(1) + to_num<Real_wp>("1.2")) > threshold)
	 || (abs(u.Value(2) - to_num<Real_wp>("0.9")) > threshold) || (abs(u.Value(3) - to_num<Real_wp>("2.7")) > threshold)
	 || (abs(u.Value(4) + to_num<Real_wp>("2.9")) > threshold) || (abs(u.Value(5) - to_num<Real_wp>("3.1")) > threshold) 
         || isnan(u.Value(0)) || isnan(u.Value(1)) || isnan(u.Value(2))
         || isnan(u.Value(3)) || isnan(u.Value(4)) || isnan(u.Value(5)))
      {
	cout << "Operator () incorrect" << endl;
	abort();
      }
    
    u.Get(12) = 1e-7;
    u.RemoveSmallEntry(1e-3);
    if ( (u.GetM() != 6) || (u.Index(0) != 1) || (u.Index(1) != 2) || (u.Index(2) != 3)
	 || (u.Index(3) != 5) || (u.Index(4) != 7) || (u.Index(5) != 8)
	 || (abs(u.Value(0) - to_num<Real_wp>("0.8")) > threshold) || (abs(u.Value(1) + to_num<Real_wp>("1.2")) > threshold)
	 || (abs(u.Value(2) - to_num<Real_wp>("0.9")) > threshold) || (abs(u.Value(3) - to_num<Real_wp>("2.7")) > threshold)
         || (abs(u.Value(4) + to_num<Real_wp>("2.9")) > threshold) || (abs(u.Value(5) - to_num<Real_wp>("3.1")) > threshold)
         || isnan(u.Value(0)) || isnan(u.Value(1)) || isnan(u.Value(2))
         || isnan(u.Value(3)) || isnan(u.Value(4)) || isnan(u.Value(5)))
      {
	cout << "RemoveSmallEntry incorrect" << endl;
	abort();
      }
    
    v.Clear();
    v.AddInteraction(4, to_num<Real_wp>("0.4"));
    Vector<int> col(3);
    Vector<Real_wp> val(3);
    col(0) = 6;
    col(1) = 0;
    col(2) = 4;
    val(0) = to_num<Real_wp>("0.8");
    val(1) = to_num<Real_wp>("2.1");
    val(2) = to_num<Real_wp>("-1.5");
    v.AddInteractionRow(3, col, val);
    
    if ( (v.GetM() != 3) || (v.Index(0) != 0) || (v.Index(1) != 4) || (v.Index(2) != 6)
	 || (abs(v.Value(0) - to_num<Real_wp>("2.1")) > threshold) || (abs(v.Value(1) + to_num<Real_wp>("1.1")) > threshold)
	 || (abs(v.Value(2) - to_num<Real_wp>("0.8")) > threshold) || isnan(v.Value(0)) || isnan(v.Value(1)) || isnan(v.Value(2)))
      {
        DISP(v.Value(0) - to_num<Real_wp>("2.1"));
        DISP(v.Value(1) + to_num<Real_wp>("1.1"));
        DISP(v.Value(2) - to_num<Real_wp>("0.8"));
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
    v.WriteText("toto.dat");
    w.ReadText("toto.dat");
    if ( (w.GetM() != 3) || (w.Index(0) != 0) || (w.Index(1) != 4) || (w.Index(2) != 6)
	 || (abs(w.Value(0) - to_num<Real_wp>("2.1")) > threshold) || (abs(w.Value(1) + to_num<Real_wp>("1.1")) > threshold)
	 || (abs(w.Value(2) - to_num<Real_wp>("0.8")) > threshold) || isnan(w.Value(0)) || isnan(w.Value(1)) || isnan(w.Value(2)) )
      {
        DISP(w); DISP(w.Value(0) - to_num<Real_wp>("2.1")); DISP(w.Value(1) + to_num<Real_wp>("1.1"));
        DISP(w.Value(2) - to_num<Real_wp>("0.8"));
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
    w.Clear();
    if (w.GetM() != 0)
      {
	cout << "Clear incorrect" << endl;
      }
    
    v.Write("totob.dat");
    w.Read("totob.dat");

    if ( (w.GetM() != 3) || (w.Index(0) != 0) || (w.Index(1) != 4) || (w.Index(2) != 6)
	 || (abs(w.Value(0) - v.Value(0)) > threshold) || (abs(w.Value(1) - v.Value(1)) > threshold)
	 || (abs(w.Value(2) - v.Value(2)) > threshold) )
      {
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
    w.Fill(2.4);
    if ( (w.GetM() != 3) || (w.Index(0) != 0) || (w.Index(1) != 4) || (w.Index(2) != 6)
	 || (abs(w.Value(0) - 2.4) > threshold) || (abs(w.Value(1) - 2.4) > threshold)
	 || (abs(w.Value(2) - 2.4) > threshold) || isnan(w.Value(0)) || isnan(w.Value(1)) || isnan(w.Value(2)))
      {
	cout << "Fill incorrect" << endl;
	abort();
      }
    
  }
  
  {
    //Real_wp threshold = 1e-14;
    // testing sparse vector of complex Real_wps
    Vector<complex<Real_wp>, VectSparse> v;
    
    v.Reallocate(4);
    if (v.GetM() != 4)
      {
	cout << "GetM incorrect" << endl;
	abort();
      }
    
    complex<Real_wp> a(1.2, 0.1), b(0.8, -0.3), c(2.7, 1.7), d(-0.3, -0.9);
    v.Index(0) = 3;
    v.Index(1) = 1;
    v.Index(2) = 5;
    v.Index(3) = 3;
    v.Value(0) = a;
    v.Value(1) = b;
    v.Value(2) = c;
    v.Value(3) = d;
    
    v.Assemble();
    if ( (v.GetM() != 3) || (v.Index(0) != 1) || (v.Index(1) != 3)
	 || (v.Index(2) != 5) || (abs(v.Value(0) - b) > threshold)
	 || (abs(v.Value(1) - a - d) > threshold) || (abs(v.Value(2) - c) > threshold) 
         || isnan(v.Value(0)) || isnan(v.Value(1)) || isnan(v.Value(2)))
      {
	cout << "Assemble incorrect" << endl;
	abort();
      }
    
    complex<Real_wp> x(3.1, 2.2), y(-2.6, 0.02);
    v.Resize(5);
    v.Index(3) = 8;
    v.Value(3) = x;
    v.Index(4) = 2;
    v.Value(4) = y;
    v.Assemble();
    
    if ( (v.GetM() != 5) || (v.Index(0) != 1) || (v.Index(1) != 2) || (v.Index(2) != 3)
	 || (v.Index(3) != 5) || (v.Index(4) != 8)
	 || (abs(v.Value(0) - b) > threshold)
	 || (abs(v.Value(1) -  y) > threshold) || (abs(v.Value(2) - a - d) > threshold)
	 || (abs(v.Value(3) - c) > threshold) || (abs(v.Value(4) - x) > threshold) 
         || isnan(v.Value(0)) || isnan(v.Value(1)) || isnan(v.Value(2)) || isnan(v.Value(3)) || isnan(v.Value(4)))
      {
	cout << "Resize/ Assemble incorrect" << endl;
	abort();
      }
    
    Vector<complex<Real_wp>, VectSparse> w;
    w.SetData(v);
    v.Nullify();
    
    if ( (w.GetM() != 5) || (w.Index(0) != 1) || (w.Index(1) != 2) || (w.Index(2) != 3)
	 || (w.Index(3) != 5) || (w.Index(4) != 8)
	 || (abs(w.Value(0) - b) > threshold)
	 || (abs(w.Value(1) - y) > threshold) || (abs(w.Value(2) - a - d) > threshold)
	 || (abs(w.Value(3) - c) > threshold) || (abs(w.Value(4) - x) > threshold) 
         || isnan(w.Value(0)) || isnan(w.Value(1)) || isnan(w.Value(2)) || isnan(w.Value(3)) || isnan(w.Value(4)))
      {
	cout << "SetData incorrect" << endl;
	abort();
      }
    
    Vector<complex<Real_wp>, VectSparse> u(w);
    
    if ( (u.GetM() != 5) || (u.Index(0) != 1) || (u.Index(1) != 2) || (u.Index(2) != 3)
	 || (u.Index(3) != 5) || (u.Index(4) != 8)
	 || (abs(u.Value(0) - b) > threshold)
	 || (abs(u.Value(1) - y) > threshold) || (abs(u.Value(2) - a - d) > threshold)
	 || (abs(u.Value(3) - c) > threshold) || (abs(u.Value(4) - x) > threshold) 
         || isnan(u.Value(0)) || isnan(u.Value(1)) || isnan(u.Value(2)) || isnan(u.Value(3)) || isnan(u.Value(4)))
      {
	cout << "Constructor incorrect" << endl;
	abort();
      }
    
    complex<Real_wp> z(-2.9, -5.3);
    complex<Real_wp> y2(1.4, 0.5);
    u.Get(2) += y2;
    u.Get(7) = z;
    
    if ( (u.GetM() != 6) || (u.Index(0) != 1) || (u.Index(1) != 2) || (u.Index(2) != 3)
	 || (u.Index(3) != 5) || (u.Index(4) != 7) || (u.Index(5) != 8)
	 || (abs(u.Value(0) - b) > threshold) || (abs(u.Value(1) - y - y2) > threshold)
	 || (abs(u.Value(2) - a-d) > threshold) || (abs(u.Value(3) - c) > threshold)
	 || (abs(u.Value(4) - z) > threshold) || (abs(u.Value(5) - x) > threshold) 
         || isnan(u.Value(0)) || isnan(u.Value(1)) || isnan(u.Value(2)) || isnan(u.Value(3)) || isnan(u.Value(4)) || isnan(u.Value(5)))
      {
	cout << "Operator () incorrect" << endl;
	abort();
      }
    
    u.Get(12) = 1e-7;
    u.RemoveSmallEntry(1e-3);
    if ( (u.GetM() != 6) || (u.Index(0) != 1) || (u.Index(1) != 2) || (u.Index(2) != 3)
	 || (u.Index(3) != 5) || (u.Index(4) != 7) || (u.Index(5) != 8)
	 || (abs(u.Value(0) - b) > threshold) || (abs(u.Value(1) - y - y2) > threshold)
	 || (abs(u.Value(2) - a - d) > threshold) || (abs(u.Value(3) - c) > threshold)
	 || (abs(u.Value(4) - z) > threshold) || (abs(u.Value(5) - x) > threshold) 
         || isnan(u.Value(0)) || isnan(u.Value(1)) || isnan(u.Value(2)) || isnan(u.Value(3)) || isnan(u.Value(4)) || isnan(u.Value(5)))
      {
	cout << "RemoveSmallEntry incorrect" << endl;
	abort();
      }
    
    v.Clear();
    v.AddInteraction(4, y);
    Vector<int> col(3);
    Vector<complex<Real_wp> > val(3);
    col(0) = 6;
    col(1) = 0;
    col(2) = 4;
    val(0) = x;
    val(1) = a;
    val(2) = b;
    v.AddInteractionRow(3, col, val);
    
    if ( (v.GetM() != 3) || (v.Index(0) != 0) || (v.Index(1) != 4) || (v.Index(2) != 6)
	 || (abs(v.Value(0) - a) > threshold) || (abs(v.Value(1) - b - y) > threshold)
	 || (abs(v.Value(2) - x) > threshold) || isnan(v.Value(0)) || isnan(v.Value(1)) || isnan(v.Value(2)))
      {
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
    v.WriteText("toto.dat");
    w.ReadText("toto.dat");
    
    if ( (w.GetM() != 3) || (w.Index(0) != 0) || (w.Index(1) != 4) || (w.Index(2) != 6)
	 || (abs(w.Value(0) - a) > threshold) || (abs(w.Value(1) - b - y) > threshold)
	 || (abs(w.Value(2) - x) > threshold) || isnan(w.Value(0)) || isnan(w.Value(1)) || isnan(w.Value(2)))
      {
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
    w.Clear();
    if (w.GetM() != 0)
      {
	cout << "Clear incorrect" << endl;
      }
    
    v.Write("totob.dat");
    w.Read("totob.dat");
    
    if ( (w.GetM() != 3) || (w.Index(0) != 0) || (w.Index(1) != 4) || (w.Index(2) != 6)
	 || (abs(w.Value(0) - a) > threshold) || (abs(w.Value(1) - b - y) > threshold)
	 || (abs(w.Value(2) - x) > threshold) )
      {
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
    w.Fill(z);
    if ( (w.GetM() != 3) || (w.Index(0) != 0) || (w.Index(1) != 4) || (w.Index(2) != 6)
	 || (abs(w.Value(0) - z) > threshold) || (abs(w.Value(1) - z) > threshold)
	 || (abs(w.Value(2) - z) > threshold) || isnan(w.Value(0)) || isnan(w.Value(1)) || isnan(w.Value(2)))
      {
	cout << "AddInteraction incorrect" << endl;
	abort();
      }
    
  }
  
  {
    // testing functions in Functions_Arrays.cxx
    int n = 30, n0 = 5, n1 = 25;
    Vector<Real_wp> x(n), x2, x3, y;
    Vector<int> permut(n), permut1(n), permut2(n);
    
    x.FillRand();
    y = x; x2 = x; x3 = x;
    permut.Fill(); permut1.Fill(); permut2.Fill();
    QuickSort(n0, n1, x);
    QuickSort(n0, n1, x2, permut);
    QuickSort(n0, n1, x3, permut1, permut2);
    for (int i = 0; i < n; i++)
      {
	if ( (x2(i) != x(i)) || (x3(i) != x(i)) || (permut1(i) != permut(i))
	     || (permut2(i) != permut(i)) || isnan(x(i)) || isnan(x2(i)) || isnan(x3(i)))
	  {
	    cout << "QuickSort incorrect" << endl;
	    abort();
	  }
	
	if ( (i >= n0) && (i < n1))
	  {
	    if ( (x(i) > x(i+1)) || (permut(i) < n0) || (permut(i) > n1)
		 || (y(permut(i)) != x(i)) || isnan(y(permut(i))))
	      {
		cout << "QuickSort incorrect" << endl;
		abort();
	      }
	  }
	else if (i == n1)
	  {
	    if ( (permut(i) < n0) || (permut(i) > n1) || (y(permut(i)) != x(i)) )
	      {
		cout << "QuickSort incorrect" << endl;
		abort();
	      }

	  }
	else
	  {
	    if ((x(i) != y(i)) || isnan(y(i)))
	      {
		cout << "QuickSort incorrect" << endl;
		abort();
	      }
	  }
      }    

    x.FillRand();
    y = x; x2 = x; x3 = x;
    permut.Fill(); permut1.Fill(); permut2.Fill();
    MergeSort(n0, n1, x);
    MergeSort(n0, n1, x2, permut);
    MergeSort(n0, n1, x3, permut1, permut2);
    for (int i = 0; i < n; i++)
      {
	if ( (x2(i) != x(i)) || (x3(i) != x(i)) || (permut1(i) != permut(i))
	     || (permut2(i) != permut(i)) || isnan(x(i)) || isnan(x2(i)) || isnan(x3(i)))
	  {
	    cout << "MergeSort incorrect" << endl;
	    abort();
	  }
	
	if ( (i >= n0) && (i < n1))
	  {
	    if ( (x(i) > x(i+1)) || (permut(i) < n0) || (permut(i) > n1)
		 || (y(permut(i)) != x(i)) || isnan(y(permut(i))))
	      {
		cout << "MergeSort incorrect" << endl;
		abort();
	      }
	  }
	else if (i == n1)
	  {
	    if ( (permut(i) < n0) || (permut(i) > n1) || (y(permut(i)) != x(i)) )
	      {
		cout << "MergeSort incorrect" << endl;
		abort();
	      }

	  }
	else
	  {
	    if (x(i) != y(i))
	      {
		cout << "MergeSort incorrect" << endl;
		abort();
	      }
	  }
      }    

    x.FillRand();
    y = x; x2 = x; x3 = x;
    permut.Fill(); permut1.Fill(); permut2.Fill();
    Sort(x);
    Sort(x2, permut);
    Sort(x3, permut1, permut2);
    for (int i = 0; i < n; i++)
      {
	if ( (x2(i) != x(i)) || (x3(i) != x(i)) || (permut1(i) != permut(i))
	     || (permut2(i) != permut(i)) || isnan(x(i)) || isnan(x2(i)) || isnan(x3(i)))
	  {
	    cout << "Sort incorrect" << endl;
	    abort();
	  }
	
	if (i < n-1)
	  {
	    if ( (x(i) > x(i+1)) || (y(permut(i)) != x(i)) || isnan(y(permut(i))))
	      {
		cout << "Sort incorrect" << endl;
		abort();
	      }
	  }
	else
	  {
	    if ( (y(permut(i)) != x(i)) )
	      {
		cout << "Sort incorrect" << endl;
		abort();
	      }

	  }
      }    
    
    int p = 15;
    Vector<int> num(n);
    Vector<bool> present(p); present.Fill(false);
    Vector<Real_wp> somme(p); somme.Fill(0);
    for (int i = 0; i < n; i++)
      {
	num(i) = rand()%p;
	x(i) = Real_wp(rand())/RAND_MAX;
	present(num(i)) = true;
	somme(num(i)) += x(i);
      }
    
    y = x;
    permut = num; permut1 = num;
    int nb = permut.GetM(); int nb2 = permut.GetM();
    Assemble(nb, permut);
    Assemble(nb2, permut1, x);
    for (int i = 0; i < nb; i++)
      {
	if ( !present(permut(i)) || (permut(i) != permut1(i)) )
	  {
	    cout << "Assemble incorrect" << endl;
	    abort();
	  }
	
	if (i < nb-1)
	  if ( (permut(i) >= permut(i+1) ) || (abs(x(i) - somme(permut(i))) > threshold) || isnan(somme(permut(i))) || isnan(x(i)))
	    {
	      cout << "Assemble incorrect" << endl;
	      abort();
	    }
      }
    
    x.Resize(nb); permut.Resize(nb);
    
    permut = num;
    permut1 = num; permut2 = num; permut2.Fill();
    RemoveDuplicate(permut);
    RemoveDuplicate(permut1, permut2);
    for (int i = 0; i < permut.GetM(); i++)
      {
	if ( !present(permut(i))  || (permut(i) != permut1(i)) || (permut1(i) != num(permut2(i)) ) )
	  {
	    cout << "RemoveDuplicate incorrect" << endl;
	    abort();
	  }

	if (i < permut.GetM()-1)
	  if (permut(i) >= permut(i+1) )
	    {
	      cout << "RemoveDuplicate incorrect" << endl;
	      abort();
	    }
      }
    

    if (!HasElement(x, x(3)))
      {
	cout << "HasElement incorrect" << endl;
	abort();
      }

    if (HasElement(x, x.GetNormInf()+Real_wp(1)))
      {
	cout << "HasElement incorrect" << endl;
	abort();
      }
  }
    
  std::remove("toto.dat");
  std::remove("totob.dat");
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
