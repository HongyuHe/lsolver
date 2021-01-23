#include "SeldonLib.hxx"

using namespace Seldon;

Real_wp threshold;

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%4;
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

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n)
{
  A.Reallocate(m, n);
  T x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        GetRandNumber(x);
        A.Set(i, j, x);
      }
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y,
		 Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - y(i)) > eps)
      return false;
  
  return true;
}

template<class T, class Prop, class Storage, class Alloc, class Storage2>
void CheckGeneralComplexEigenvalue(Matrix<complex<T>, Prop, Storage, Alloc>& A,
				   int n, const Storage2& stor)
{
  GenerateRandomMatrix(A, n, n);
  
  // checking GetEigenvaluesEigenvectors
  Matrix<complex<T>, Prop, Storage, Alloc> B(A), C, D;
  Vector<complex<T> > lambda, x(n), y(n);
  Matrix<complex<T>, General, Storage2> V, V2;
  GetEigenvaluesEigenvectors(A, lambda, V);
  
  for (int j = 0; j < n; j++)
    {
      // checking if A x = lambda x
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);
      
      // y = A x
      Mlt(B, x, y);

      Real_wp err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(y(i) - lambda(j)*x(i)), 2.0);
      
      err = sqrt(err);
      if (abs(lambda(j)) > threshold)
        if (err > threshold*abs(lambda(j)))
          {
            cout << "Error on eigenvalue " << lambda(j) << endl;
            cout << "Error = " << err << endl;
            abort();
          }
    }

  // checking GetEigenvalues
  Vector<complex<T> > w;
  A = B;
  GetEigenvalues(A, w);
  
  for (int i = 0; i < n; i++)
    if (abs(lambda(i) - w(i)) > threshold)
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      

  // checking generalized eigenvalues
  GenerateRandomMatrix(A, n, n);
  GenerateRandomMatrix(B, n, n);
  Vector<complex<T> > Mx(n), alpha, beta; Mx.Fill(0);
  C = A; D = B;
  GetEigenvaluesEigenvectors(A, B, alpha, beta, V);
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);      
      
      // y = A x, Mx = B x
      Mlt(C, x, y);
      Mlt(D, x, Mx);

      T err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(beta(j)*y(i) - alpha(j)*Mx(i)), 2.0);
      
      err = sqrt(err);
      if (err > threshold)
	{
	  cout << "Error on eigenvalue " << alpha(j)/beta(j) << endl;
	  cout << "Error = " << err << endl;
	  abort();
	}
    }
  
  A = C; B = D;
  Vector<complex<T> > wdiv;
  GetEigenvalues(A, B, w, wdiv);
  for (int i = 0; i < n; i++)
    if ( (abs(alpha(i) - w(i)) > threshold) || (abs(beta(i) - wdiv(i)) > threshold))
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      

}


template<class T, class Prop, class Storage, class Alloc, class Storage2>
void CheckSymmetricComplexEigenvalue(Matrix<complex<T>, Prop, Storage, Alloc>& A,
				     int n, const Storage2& stor)
{
  GenerateRandomMatrix(A, n, n);
  
  // checking GetEigenvaluesEigenvectors
  Matrix<complex<T>, Prop, Storage, Alloc> B(A), C, D;
  Vector<complex<T> > lambda, x(n), y(n);
  Matrix<complex<T>, General, Storage2> V, V2;
  GetEigenvaluesEigenvectors(A, lambda, V);
  
  for (int j = 0; j < n; j++)
    {
      // checking if A x = lambda x
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);
      
      // y = A x
      Mlt(B, x, y);

      Real_wp err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(y(i) - lambda(j)*x(i)), 2.0);
      
      err = sqrt(err);
      if (abs(lambda(j)) > threshold)
        if (err > threshold*abs(lambda(j)))
          {
            cout << "Error on eigenvalue " << lambda(j) << endl;
            cout << "Error = " << err << endl;
            abort();
          }
    }

  // checking GetEigenvalues
  Vector<complex<T> > w;
  A = B;
  GetEigenvalues(A, w);
  
  for (int i = 0; i < n; i++)
    if (abs(lambda(i) - w(i)) > threshold)
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      

  // checking generalized eigenvalues
  GenerateRandomMatrix(A, n, n);
  GenerateRandomMatrix(B, n, n);
  Vector<complex<T> > Mx(n), alpha, beta; Mx.Fill(0);
  C = A; D = B;
  GetEigenvaluesEigenvectors(A, B, alpha, beta, V);
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);      
      
      // y = A x, Mx = B x
      Mlt(C, x, y);
      Mlt(D, x, Mx);

      T err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(beta(j)*y(i) - alpha(j)*Mx(i)), 2.0);
      
      err = sqrt(err);
      if (err > threshold)
	{
	  cout << "Error on eigenvalue " << alpha(j)/beta(j) << endl;
	  cout << "Error = " << err << endl;
	  abort();
	}
    }
  
  A = C; B = D;
  Vector<complex<T> > wdiv;
  GetEigenvalues(A, B, w, wdiv);
  for (int i = 0; i < n; i++)
    if ( (abs(alpha(i) - w(i)) > threshold) || (abs(beta(i) - wdiv(i)) > threshold))
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }

}


template<class T, class Prop, class Storage, class Alloc>
void CheckGeneralRealEigenvalue(Matrix<T, Prop, Storage, Alloc>& A, int n)
{
  GenerateRandomMatrix(A, n, n);
  
  // checking GetEigenvaluesEigenvectors
  Matrix<T, Prop, Storage, Alloc> B(A), C, D;
  Vector<T> lambda_r, lambda_i;
  Matrix<T, General, Storage> V, V2;
  GetEigenvaluesEigenvectors(A, lambda_r, lambda_i, V);
  Vector<T> X(n), Xi(n), Y(n), Yi(n);
  X.Fill(0); Xi.Fill(0);
  Y.Fill(0); Yi.Fill(0);
  
  int i = 0;
  while (i < lambda_r.GetM())
    {
      // on regarde si la valeur propre est reelle ou si elle fait
      // partie d'une paire de valeurs propres conjuguees
      bool eigen_pair = false;
      if (i < lambda_r.GetM()-1)
        {
          if ( (lambda_r(i) == lambda_r(i+1)) && (lambda_i(i) == -lambda_i(i+1)))
            eigen_pair = true;
        }
      
      T err = 0;
      if (eigen_pair)
        {
	  // paire, X + i Xi et X - i Xi sont les vecteurs propres
	  // associes a (lambda_r, lambda_i) et (lambda_r, -lambda_i)
          for (int j = 0; j < n; j++)
            {
              X(j) = V(j, i);
              Xi(j) = V(j, i+1);
            }
          
          Mlt(B, X, Y);
          Mlt(B, Xi, Yi);
          for (int j = 0; j < n; j++)
            err += pow(abs(complex<T>(Y(j), Yi(j))
                           - complex<T>(lambda_r(i), lambda_i(i))
                           *complex<T>(X(j), Xi(j)) ), 2.0);
          
          err = sqrt(err);
        }
      else
        {
	  // valeur propre reelle
          for (int j = 0; j < n; j++)
            X(j) = V(j, i);
          
          Mlt(B, X, Y);
          for (int j = 0; j < n; j++)
            err += pow(Y(j) - lambda_r(i)*X(j), 2.0);
          
          err = sqrt(err);
        }
      
      if (err > threshold*abs(complex<T>(lambda_r(i), lambda_i(i))))
        {
          cout << "Error on eigenvalue " << lambda_r(i) << endl;
          cout << "Error = " << err << endl;
          abort();
        }
      
      if (eigen_pair)
        i += 2;
      else
        i++;
    }
  
  Vector<T> wr, wi;
  A = B;
  GetEigenvalues(A, wr, wi);
  for (int i = 0; i < n; i++)
    if ( (abs(lambda_r(i) - wr(i)) > threshold) ||
	 (abs(lambda_i(i) - wi(i)) > threshold) )
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      
 
  // checking generalized eigenvalues
  GenerateRandomMatrix(A, n, n);
  GenerateRandomMatrix(B, n, n);
  Vector<T> Mx(n), Mxi(n), alpha_r, alpha_i, beta;
  Mx.Fill(0); Mxi.Fill(0);
  C = A; D = B;
  GetEigenvaluesEigenvectors(A, B, alpha_r, alpha_i, beta, V);
  while (i < alpha_r.GetM())
    {
      // on regarde si la valeur propre est reelle ou si elle fait
      // partie d'une paire de valeurs propres conjuguees
      bool eigen_pair = false;
      if (i < lambda_r.GetM()-1)
        {
          if ( (lambda_r(i) == lambda_r(i+1)) && (lambda_i(i) == -lambda_i(i+1)))
            eigen_pair = true;
        }
      
      T err = 0;
      if (eigen_pair)
        {
	  // paire, X + i Xi et X - i Xi sont les vecteurs propres
	  // associes a (lambda_r, lambda_i) et (lambda_r, -lambda_i)
          for (int j = 0; j < n; j++)
            {
              X(j) = V(j, i);
              Xi(j) = V(j, i+1);
            }
          
          Mlt(C, X, Y);
          Mlt(C, Xi, Yi);
          Mlt(D, X, Mx);
          Mlt(D, Xi, Mxi);
          for (int j = 0; j < n; j++)
            err += pow(abs(beta(i)*complex<T>(Y(j), Yi(j))
                           - complex<T>(alpha_r(i), alpha_i(i))
                           *complex<T>(Mx(j), Mxi(j)) ), 2.0);
          
          err = sqrt(err);
        }
      else
        {
	  // valeur propre reelle
          for (int j = 0; j < n; j++)
            X(j) = V(j, i);
          
          Mlt(C, X, Y);
          Mlt(D, X, Mx);
          for (int j = 0; j < n; j++)
            err += pow(beta(i)*Y(j) - lambda_r(i)*Mx(j), 2.0);
          
          err = sqrt(err);
        }
      
      if (err > threshold)
        {
          cout << "Error on eigenvalue " << alpha_r(i) << endl;
          cout << "Error = " << err << endl;
          abort();
        }
      
      if (eigen_pair)
        i += 2;
      else
        i++;
    }
  
  A = C; B = D;
  Vector<T> wdiv;
  GetEigenvalues(A, B, lambda_r, lambda_i, wdiv);
  for (int i = 0; i < n; i++)
    if ( (abs(alpha_r(i) - lambda_r(i)) > threshold) ||
	 (abs(alpha_i(i) - lambda_i(i)) > threshold) ||
	 (abs(beta(i) - wdiv(i)) > threshold))
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }        

}

template<class T, class Prop, class Storage, class Alloc, class Storage2>
void CheckSymmetricRealEigenvalue(Matrix<T, Prop, Storage, Alloc>& A, int n,
				  const Storage2& stor)
{
  GenerateRandomMatrix(A, n, n);
  
  // checking GetEigenvaluesEigenvectors
  Matrix<T, Prop, Storage, Alloc> B(A);
  Vector<T> lambda, x(n), y(n);
  Matrix<T, General, Storage2> V, V2;
  GetEigenvaluesEigenvectors(A, lambda, V);
  
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);      
      
      // y = A x
      Mlt(B, x, y);

      T err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(y(i) - lambda(j)*x(i)), 2.0);
      
      err = sqrt(err);
      if (abs(lambda(j)) > threshold)
        if (err > threshold*abs(lambda(j)))
          {
            cout << "Error on eigenvalue " << lambda(j) << endl;
            cout << "Error = " << err << endl;
            abort();
          }
    }

  Vector<T> w;
  A = B;
  GetEigenvalues(A, w);
  for (int i = 0; i < n; i++)
    if (abs(lambda(i) - w(i)) > threshold)
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      
  
  Matrix<T, Prop, Storage, Alloc> C, D;
  
  // checking generalized eigenvalues
  GenerateRandomMatrix(A, n, n);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      {
        B(i, j) = 0;
        for (int k = 0; k < n; k++)
          B(i, j) += A(i, k)*A(j, k);
      }

  Vector<T> Mx(n); Mx.Fill(0);
  C = A; D = B;
  GetEigenvaluesEigenvectors(A, B, lambda, V);
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);      
      
      // y = A x, Mx = B x
      Mlt(C, x, y);
      Mlt(D, x, Mx);

      T err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(y(i) - lambda(j)*Mx(i)), 2.0);
      
      err = sqrt(err);
      if (abs(lambda(j)) > threshold)
        if (err > threshold*abs(lambda(j)))
          {
            cout << "Error on eigenvalue " << lambda(j) << endl;
            cout << "Error = " << err << endl;
            abort();
          }
    }
  
  A = C; B = D;
  GetEigenvalues(A, B, w);
  for (int i = 0; i < n; i++)
    if (abs(lambda(i) - w(i)) > threshold)
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      

}

template<class T, class Prop, class Storage, class Alloc, class Storage2>
void CheckHermitianEigenvalue(Matrix<complex<T>, Prop, Storage, Alloc>& A, int n,
			      const Storage2& stor)
{
  GenerateRandomMatrix(A, n, n);
  for (int i = 0; i < n; i++)
    A.Get(i, i) = real(A(i, i));
  
  // checking GetEigenvaluesEigenvectors
  Matrix<complex<T>, Prop, Storage, Alloc> B(A), C, D;
  Vector<T> lambda;
  Vector<complex<T> > x(n), y(n);
  Matrix<complex<T>, General, Storage2> V;
  GetEigenvaluesEigenvectors(A, lambda, V);
  
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);      
      
      // y = A x
      Mlt(B, x, y);

      T err = 0;
      for (int i = 0; i < n; i++)
        err += absSquare(y(i) - lambda(j)*x(i));
      
      err = sqrt(err);
      if (abs(lambda(j)) > threshold)
        if (err > threshold*abs(lambda(j)))
          {
            cout << "Error on eigenvalue " << lambda(j) << endl;
            cout << "Error = " << err << endl;
            abort();
          }
    }

  Vector<T> w;
  A = B;
  GetEigenvalues(A, w); 
  for (int i = 0; i < n; i++)
    if (abs(lambda(i) - w(i)) > threshold)
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      
  
  // checking generalized eigenvalues
  Matrix<complex<T> > L, LLc;
  GenerateRandomMatrix(A, n, n);
  for (int i = 0; i < n; i++)
    A.Get(i, i) = real(A(i, i));

  GenerateRandomMatrix(L, n, n);
  LLc.Reallocate(n, n);
  MltAdd(complex<T>(1.0, 0.0), SeldonNoTrans, L, SeldonConjTrans, L, complex<T>(0.0, 0.0), LLc);  
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      B.Get(i, j) = LLc(i, j);
  
  Vector<complex<T> > Mx(n); Mx.Fill(0);
  C = A; D = B;
  GetEigenvaluesEigenvectors(A, B, lambda, V);
  for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
        x(i) = V(i, j);      
      
      // y = A x, Mx = B x
      Mlt(C, x, y);
      Mlt(D, x, Mx);

      T err = 0;
      for (int i = 0; i < n; i++)
        err += pow(abs(y(i) - lambda(j)*Mx(i)), 2.0);
      
      err = sqrt(err);
      if (err > 10.0*Norm2(y)*threshold)
	{
	  cout << "Error on eigenvalue " << lambda(j) << endl;
	  cout << "Error = " << err << endl;
	  abort();
	}
    }
  
  A = C; B = D;
  GetEigenvalues(A, B, w);
  for (int i = 0; i < n; i++)
    if (abs(lambda(i) - w(i)) > threshold)
      {
        cout << "GetEigenvalues incorrect" << endl;
        abort();
      }      

}

template<class T, class Prop, class Storage, class Alloc>
bool CheckIdentity(const Matrix<T, Prop, Storage, Alloc>& A)
{
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == j)
          {
            if (abs(A(i, j) - 1.0) > threshold)
              return false;
          }
        else
          {
            if (abs(A(i, j)) > threshold)
              return false;
          }
      }
  
  return true;
}

template<class T1, class Prop1, class Storage1, class Allocator1,
         class T2, class Prop2, class Storage2, class Allocator2>
bool EqualMatrix(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
                 const Matrix<T2, Prop2, Storage2, Allocator2>& B,
		 Real_wp eps = threshold)
{
  if ( (A.GetM() != B.GetM())  || (A.GetN() != B.GetN()) )
    {
      DISP(A.GetM()); DISP(A.GetN());
      DISP(B.GetM()); DISP(B.GetN());
      return false;
    }
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if (abs(A(i, j) - B(i, j)) > eps)
        {
	  DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
	  return false;
	}
  
  return true;
}


template<class T, class Prop, class Storage, class Allocator>
void CheckSylvester(Matrix<T, Prop, Storage, Allocator>& A, int n)
{
  Matrix<T, Prop, Storage, Allocator> B, Q, Z, AA, BB, Id, C, D, E, X, X0;
  GenerateRandomMatrix(A, n, n);
  GenerateRandomMatrix(B, n, n);

  T zero, one;
  SetComplexZero(zero);
  SetComplexOne(one);
    
  // checking GetHessenberg
  AA = A; BB = B;
  GetHessenberg(AA, Q);
  Id.Reallocate(n, n);
  MltAdd(one, SeldonNoTrans, Q, SeldonConjTrans, Q, zero, Id);
  
  if (!CheckIdentity(Id))
    {
      cout << "GetHessenberg incorrect" << endl;
      abort();
    }

  C.Reallocate(n, n);
  C.Fill(zero);
  Mlt(A, Q, C);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C, zero, Id);
  if (!EqualMatrix(Id, AA))
    {
      cout << "GetHessenberg incorrect" << endl;
      abort();
    }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < i-1; j++)
      {
	if ( abs(AA(i, j)) > threshold )
	  {
	    cout << "GetHessenberg incorrect" << endl;
	    abort();
	  }
      }
    
  AA = A;
  GetHessenberg(AA, BB, Q, Z);
  
  Id.Reallocate(n, n);
  MltAdd(one, SeldonNoTrans, Q, SeldonConjTrans, Q, zero, Id);
  
  if (!CheckIdentity(Id))
    {
      cout << "GetHessenberg incorrect" << endl;
      abort();
    }

  MltAdd(one, SeldonNoTrans, Z, SeldonConjTrans, Z, zero, Id);
  
  if (!CheckIdentity(Id))
    {
      cout << "GetHessenberg incorrect" << endl;
      abort();
    }
  
  C.Reallocate(n, n);
  C.Fill(zero);
  Mlt(A, Z, C);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C, zero, Id);
  if (!EqualMatrix(Id, AA))
    {
      cout << "GetHessenberg incorrect" << endl;
      abort();
    }

  Mlt(B, Z, C);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C, zero, Id);
  if (!EqualMatrix(Id, BB))
    {
      cout << "GetHessenberg incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      {
	if ( (j < i-1) && (abs(AA(i, j)) > threshold) )
	  {
	    cout << "GetHessenberg incorrect" << endl;
	    abort();
	  }
	
	if (abs(BB(i, j)) > threshold )
          {
	    cout << "GetHessenberg incorrect" << endl;
	    abort();
	  }
      }
  
  // checking GetQZ
  AA = A; BB = B;
  GetQZ(AA, BB, Q, Z);
  
  Id.Reallocate(n, n);
  MltAdd(one, SeldonNoTrans, Q, SeldonConjTrans, Q, zero, Id);
  
  if (!CheckIdentity(Id))
    {
      cout << "GetQZ incorrect" << endl;
      abort();
    }

  MltAdd(one, SeldonNoTrans, Z, SeldonConjTrans, Z, zero, Id);
  
  if (!CheckIdentity(Id))
    {
      cout << "GetQZ incorrect" << endl;
      abort();
    }
  
  C.Reallocate(n, n);
  C.Fill(zero);
  Mlt(A, Z, C);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C, zero, Id);
  if (!EqualMatrix(Id, AA))
    {
      cout << "GetQZ incorrect" << endl;
      abort();
    }

  Mlt(B, Z, C);
  MltAdd(one, SeldonConjTrans, Q, SeldonNoTrans, C, zero, Id);
  if (!EqualMatrix(Id, BB))
    {
      cout << "GetQZ incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      {
	if ((abs(AA(i, j)) > threshold) || (abs(BB(i, j)) > threshold) )
	  {
	    bool block22 = false;
	    if (j == i-1)
	      {
		block22 = true;
		if (j > 0)
		  if (abs(AA(i-1, j-1)) > threshold)
		    block22 = false;
		
		if (i < n-1)
		  if (abs(AA(i+1, j+1)) > threshold)
		    block22 = false;		  
	      }
	    
	    if (!block22)
	      {
		cout << "GetQZ incorrect" << endl;
		abort();
	      }
	  }
      }
  
  // checking SolveHessenberg, SolveHessenbergTwo
  GenerateRandomMatrix(A, n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i-1; j++)
      A(i, j) = zero;
  
  Vector<T> x(n), b(n), y;
  GenerateRandomVector(x, n);
  Mlt(A, x, b);
  
  y = x;
  x = b;
  SolveHessenberg(A, x);
  
  if (!EqualVector(x, y))
    {
      cout << "SolveHessenberg incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(A, n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i-2; j++)
      A(i, j) = zero;
  
  GenerateRandomVector(x, n);
  Mlt(A, x, b);
  
  y = x;
  x = b;
  SolveHessenbergTwo(A, x);
  
  if (!EqualVector(x, y))
    {
      cout << "SolveHessenbergTwo incorrect" << endl;
      abort();
    }
  
  // checking SolveSylvester
  GenerateRandomMatrix(A, n, n);
  GenerateRandomMatrix(B, n, n);
  GenerateRandomMatrix(C, n, n);
  GenerateRandomMatrix(D, n, n);
  GenerateRandomMatrix(X, n, n);
  X0 = X;
  
  E.Reallocate(n, n); E.Fill(zero);
  Mlt(A, X, Q);
  MltAdd(one, SeldonNoTrans, Q, SeldonConjTrans, B, zero, E);
  Mlt(C, X, Q);
  MltAdd(one, SeldonNoTrans, Q, SeldonConjTrans, D, one, E);
  
  SolveSylvester(A, B, C, D, E);
  if (!EqualMatrix(X, E, 10.0*threshold))
    {
      cout << "SolveSylvester incorrect" << endl;
      abort();
    }
  
  // checking GetSVD, GetPseudoInverse
  int m = n+2;
  GenerateRandomMatrix(A, m, n);
  Vector<Real_wp> sigma;
  Matrix<T, Prop, Storage, Allocator> U, V;
  D = A;
  GetSVD(A, sigma, U, V);
  
  Id.Reallocate(m, m);
  Id.Fill(zero);
  MltAdd(one, SeldonNoTrans, U, SeldonConjTrans, U, zero, Id);
  if (!CheckIdentity(Id))
    {
      cout << "GetSVD" << endl;
      abort();
    }

  Id.Reallocate(n, n);
  Id.Fill(zero);
  MltAdd(one, SeldonNoTrans, V, SeldonConjTrans, V, zero, Id);
  if (!CheckIdentity(Id))
    {
      cout << "GetSVD" << endl;
      abort();
    }
  
  B.Reallocate(m, n);
  B.Fill(zero);
  for (int k = 0; k < min(m, n); k++)
    for (int j = 0; j < m; j++)
      U(j, k) *= sigma(k);
  
  V.Resize(m, n);
  for (int i = 0; i < n; i++)
    for (int j = n; j < m; j++)
      V(j, i) = 0.0;
  
  MltAdd(one, SeldonNoTrans, U, SeldonNoTrans, V, zero, B);
  if (!EqualMatrix(D, B))
    {
      cout << "GetSVD incorrect" << endl;
      abort();
    }
  
  A = D;
  GetPseudoInverse(D, threshold);
  Id.Reallocate(n, n);
  Mlt(D, A, Id);
  if (!CheckIdentity(Id))
    {
      cout << "GetPseudoInverse" << endl;
      abort();
    }
  
  m = n-2;
  GenerateRandomMatrix(A, m, n);
  D = A;
  GetSVD(A, sigma, U, V);
  
  Id.Reallocate(m, m);
  Id.Fill(zero);
  MltAdd(one, SeldonNoTrans, U, SeldonConjTrans, U, zero, Id);
  if (!CheckIdentity(Id))
    {
      cout << "GetSVD" << endl;
      abort();
    }

  Id.Reallocate(n, n);
  Id.Fill(zero);
  MltAdd(one, SeldonNoTrans, V, SeldonConjTrans, V, zero, Id);
  if (!CheckIdentity(Id))
    {
      cout << "GetSVD" << endl;
      abort();
    }
  
  B.Reallocate(m, n);
  B.Fill(zero);
  for (int k = 0; k < min(m, n); k++)
    for (int j = 0; j < m; j++)
      U(j, k) *= sigma(k);
  
  V.Resize(m, n);
  for (int i = 0; i < n; i++)
    for (int j = n; j < m; j++)
      V(j, i) = 0.0;
  
  MltAdd(one, SeldonNoTrans, U, SeldonNoTrans, V, zero, B);
  if (!EqualMatrix(D, B))
    {
      cout << "GetSVD incorrect" << endl;
      abort();
    }

  A = D;
  GetPseudoInverse(D, threshold);
  Id.Reallocate(m, m);
  Mlt(A, D, Id);
  if (!CheckIdentity(Id))
    {
      cout << "GetPseudoInverse" << endl;
      abort();
    }
  
}

int main(int argc, char** argv)
{
  bool all_test = true;
  
  threshold = 1e-11;
  
  RowMajor row_major;
  ColMajor col_major;
  
  {
    Matrix<Complex_wp, General, RowMajor> A;
    int n = 20;
    CheckGeneralComplexEigenvalue(A, n, row_major);
  }
  
  {
    Matrix<Complex_wp, General, ColMajor> A;
    int n = 20;
    CheckGeneralComplexEigenvalue(A, n, col_major);
  }

  {
    Matrix<Real_wp, General, RowMajor> A;
    int n = 20;
    CheckGeneralRealEigenvalue(A, n);
  }

  {
    Matrix<Real_wp, General, ColMajor> A;
    int n = 20;
    CheckGeneralRealEigenvalue(A, n);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY  
  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    int n = 7;
    CheckSymmetricRealEigenvalue(A, n, col_major);
  }    
    
  {
    Matrix<Complex_wp, Symmetric, ColSymPacked> A;
    int n = 14;
    CheckSymmetricComplexEigenvalue(A, n, col_major);
  }    

  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    int n = 7;
    CheckSymmetricRealEigenvalue(A, n, col_major);
  }    

  {
    Matrix<Complex_wp, Symmetric, ColSym> A;
    int n = 14;
    CheckSymmetricComplexEigenvalue(A, n, col_major);
  }    
#endif
  
  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    int n = 5;
    CheckSymmetricRealEigenvalue(A, n, row_major);
  }    

  {
    Matrix<Complex_wp, Symmetric, RowSymPacked> A;
    int n = 11;
    CheckSymmetricComplexEigenvalue(A, n, row_major);
  }    

  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    int n = 5;
    CheckSymmetricRealEigenvalue(A, n, row_major);
  }    

  {
    Matrix<Complex_wp, Symmetric, RowSym> A;
    int n = 11;
    CheckSymmetricComplexEigenvalue(A, n, row_major);
  }    

  {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    int n = 12;
    CheckHermitianEigenvalue(A, n, row_major);
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    int n = 14;
    CheckHermitianEigenvalue(A, n, col_major);
  }    
#endif

  {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    int n = 12;
    CheckHermitianEigenvalue(A, n, row_major);
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    int n = 14;
    CheckHermitianEigenvalue(A, n, col_major);
  }    
#endif
  
  {
    Matrix<Complex_wp, General, ColMajor> A;
    int n = 13;
    CheckSylvester(A, n);
  }    

  {
    Matrix<Real_wp, General, ColMajor> A;
    int n = 13;
    CheckSylvester(A, n);
  }    

  {
    Matrix<Complex_wp, General, RowMajor> A;
    int n = 9;
    CheckSylvester(A, n);
  }    

  {
    Matrix<Real_wp, General, RowMajor> A;
    int n = 11;
    CheckSylvester(A, n);
  }    


  if (all_test)
    cout << "All tests passed successfully" << endl;

  return 0;
}
