#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<double> VectReal_wp;

Real_wp threshold = 1e-12;

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  int type = rand()%3;
  //int type = 2;
  if (type == 0)
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void FillRand(Vector<T> & x)
{
  for (int i = 0; i < x.GetM(); i++)
    GetRand(x(i));  
}

template<class T, class Prop1, class Storage1, class T2, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop1, Storage1>& A,
		 const Matrix<T2, Prop2, Storage2>& B)
{
  if ( (A.GetM() != B.GetM()) || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ( (abs(A(i, j) - B(i, j) ) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          //abort();
          return false;
        }
  
  return true;
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, RowComplexSparse>& A,
                    Matrix<T2, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ColComplexSparse>& A,
                    Matrix<T2, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<complex<T>, Prop, ArrayRowComplexSparse>& A,
                    Matrix<T2, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w, z2, w2;
  SetComplexZero(zero);
  complex<T> czero(0, 0);
  GetRand(x); GetRand(y);
  GetRand(z); GetRand(w);
  GetRand(z2); GetRand(w2);

  int n = Ad.GetN();
  for (int i = 0; i < Ad.GetM(); i++)
    {
      int nb = 0, nb_i = 0;
      for (int j = 0; j < n; j++)
        {
          if (real(Ad(i, j)) != zero)
            nb++;
          
          if (imag(Ad(i, j)) != zero)
            nb_i++;
        }
      
      if (nb != A.GetRealRowSize(i))
        {
          cout << "GetRealRowSize incorrect" << endl;
          abort();
        }

      if (nb_i != A.GetImagRowSize(i))
        {
          cout << "GetImagRowSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb; j++)
        {
          int jcol = A.IndexReal(i, j);
          if ((real(Ad(i, jcol)) != A.ValueReal(i, j)) || isnan(real(Ad(i, jcol))) || isnan(A.ValueReal(i, j)))
            {
              cout << "IndexReal/ValueReal incorrect" << endl;
              abort();
            }
        }

      for (int j = 0; j < nb_i; j++)
        {
          int jcol = A.IndexImag(i, j);
          if ((imag(Ad(i, jcol)) != A.ValueImag(i, j)) || isnan(imag(Ad(i, jcol))) || isnan(A.ValueImag(i, j)))
            {
              cout << "IndexImag/ValueImag incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateRealRow(2, 2);
  A.IndexReal(2, 0) = 2; A.IndexReal(2, 1) = 5; 
  A.ValueReal(2, 0) = x; A.ValueReal(2, 1) = y; 
  A.ReallocateImagRow(2, 1);
  A.IndexImag(2, 0) = 3;
  A.ValueImag(2, 0) = w;
  Ad(2,0) = czero; Ad(2,1) = czero; Ad(2,2) = complex<T>(x, zero);
  Ad(2,3) = complex<T>(zero, w); Ad(2,4) = czero;
  Ad(2,5) = complex<T>(y, zero); Ad(2,6) = czero;
  
  A.ClearRealRow(3);
  for (int i = 0; i < n; i++)
    Ad(3, i) = complex<T>(zero, imag(Ad(3, i)));

  A.ClearImagRow(1);
  for (int i = 0; i < n; i++)
    Ad(1, i) = complex<T>(real(Ad(1, i)), 0);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateRealRow/ClearRealRow incorrect" << endl;
      abort();
    }
  
  A.ResizeRealRow(2, 4);
  A.IndexReal(2, 2) = 1; A.IndexReal(2, 3) = 4;
  A.ValueReal(2, 2) = z; A.ValueReal(2, 3) = w;
  A.AssembleRealRow(2);

  A.ResizeImagRow(1, 3);
  A.IndexImag(1, 0) = 0; A.IndexImag(1, 1) = 2; A.IndexImag(1, 2) = 6;
  A.ValueImag(1, 0) = x; A.ValueImag(1, 1) = z2; A.ValueImag(1, 2) = w2;
  A.AssembleImagRow(1);

  Ad(2, 1) = complex<T>(z, imag(Ad(2, 1))); Ad(2, 4) = complex<T>(w, imag(Ad(2, 4)));
  Ad(1, 0) = complex<T>(real(Ad(1, 0)), x); Ad(1, 2) = complex<T>(real(Ad(1, 2)), z2); Ad(1, 6) = complex<T>(real(Ad(1, 6)), w2);
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeRow incorrect" << endl;
      abort();
    }
  
  A.SwapRealRow(1, 3);
  A.SwapImagRow(2, 0);
  for (int j = 0; j < n; j++)
    {
      x = real(Ad(1, j));
      Ad(1, j) = complex<T>(real(Ad(3, j)), imag(Ad(1, j)));
      Ad(3, j) = complex<T>(x, imag(Ad(3, j)));

      x = imag(Ad(2, j));
      Ad(2, j) = complex<T>(real(Ad(2, j)), imag(Ad(0, j)));
      Ad(0, j) = complex<T>(real(Ad(0, j)), x);
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "SwapRealRow/SwapImagRow incorrect" << endl;
      abort();
    }
  
  IVect new_ind(3);
  x = A.ValReal(0, 0); y = A.ValReal(0, 2);
  new_ind(0) = 1; new_ind(1) = 4; new_ind(2) = 5;
  A.ReplaceRealIndexRow(0, new_ind);
  Ad(0, 0) = complex<T>(zero, imag(Ad(0, 0))); Ad(0, 2) = complex<T>(zero, imag(Ad(0, 2)));
  Ad(0, 1) = complex<T>(x, imag(Ad(0, 1))); Ad(0, 4) = complex<T>(y, imag(Ad(0, 4)));
  
  new_ind(0) = 0; new_ind(1) = 2; new_ind(2) = 3;
  y = A.ValImag(1, 6);
  A.ReplaceImagIndexRow(1, new_ind);
  Ad(1, 6) = complex<T>(real(Ad(1, 6)), zero); Ad(1, 3) = complex<T>(real(Ad(1, 3)), y);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReplaceIndexRealRow incorrect" << endl;
      abort();
    }

}


template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<complex<T>, Prop, ArrayColComplexSparse>& A,
                    Matrix<T2, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w, x2, y2;
  SetComplexZero(zero);
  complex<T> czero(0, 0);
  GetRand(x); GetRand(y);
  GetRand(x2); GetRand(y2);
  GetRand(z); GetRand(w);

  int m = Ad.GetM();
  for (int i = 0; i < Ad.GetN(); i++)
    {
      int nb_r = 0, nb_i = 0;
      for (int j = 0; j < m; j++)
        {
          if (real(Ad(j, i)) != zero)
            nb_r++;

          if (imag(Ad(j, i)) != zero)
            nb_i++;
        }
      
      if (nb_r != A.GetRealColumnSize(i))
        {
          cout << "GetRealColumnSize incorrect" << endl;
          abort();
        }

      if (nb_i != A.GetImagColumnSize(i))
        {
          cout << "GetImagColumnSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb_r; j++)
        {
          int jrow = A.IndexReal(i, j);
          if ((real(Ad(jrow, i)) != A.ValueReal(i, j)) || isnan(real(Ad(jrow, i))) || isnan(A.ValueReal(i, j)))
            {
              cout << "IndexReal/ValueReal incorrect" << endl;
              abort();
            }
        }
      
      for (int j = 0; j < nb_i; j++)
        {
          int jrow = A.IndexImag(i, j);
          if ((imag(Ad(jrow, i)) != A.ValueImag(i, j)) || isnan(imag(Ad(jrow, i))) || isnan(A.ValueImag(i, j)))
            {
              cout << "IndexImag/ValueImag incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateRealColumn(2, 2);
  A.IndexReal(2, 0) = 2; A.IndexReal(2, 1) = 3; 
  A.ValueReal(2, 0) = x; A.ValueReal(2, 1) = y; 
  Ad(0, 2) = complex<T>(zero, imag(Ad(0, 2))); Ad(1, 2) = complex<T>(zero, imag(Ad(1, 2)));
  Ad(2, 2) = complex<T>(x, imag(Ad(2, 2))); Ad(3, 2) = complex<T>(y, imag(Ad(3,2)));
  
  A.ClearRealColumn(5);
  for (int i = 0; i < m; i++)
    Ad(i, 5) = complex<T>(zero, imag(Ad(i, 5)));

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateRealColumn/ClearRealColumn incorrect" << endl;
      abort();
    }

  A.ReallocateImagColumn(4, 2);
  A.IndexImag(4, 0) = 0; A.IndexImag(4, 1) = 3; 
  A.ValueImag(4, 0) = x2; A.ValueImag(4, 1) = y2; 
  Ad(0, 4) = complex<T>(real(Ad(0, 4)), x2); Ad(1, 4) = complex<T>(real(Ad(1, 4)), zero);
  Ad(2, 4) = complex<T>(real(Ad(2, 4)), zero);
  Ad(3, 4) = complex<T>(real(Ad(3, 4)), y2);
  
  A.ClearImagColumn(6);
  for (int i = 0; i < m; i++)
    Ad(i, 6) = complex<T>(real(Ad(i, 6)), zero);

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateImagColumn/ClearImagColumn incorrect" << endl;
      abort();
    }
  
  A.ResizeRealColumn(2, 3);
  A.IndexReal(2, 2) = 0;
  A.ValueReal(2, 2) = z;
  A.AssembleRealColumn(2);
  Ad(0, 2) = complex<T>(z, imag(Ad(0, 2)));
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeRealColumn incorrect" << endl;
      abort();
    }

  A.ResizeImagColumn(4, 3);
  A.IndexImag(4, 2) = 1;
  A.ValueImag(4, 2) = w;
  A.AssembleImagColumn(4);
  Ad(1, 4) = complex<T>(real(Ad(1, 4)), w);
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeImagColumn incorrect" << endl;
      abort();
    }
  
  A.SwapRealColumn(2, 5);
  for (int j = 0; j < m; j++)
    {
      x = real(Ad(j, 2));
      Ad(j, 2) = complex<T>(real(Ad(j, 5)), imag(Ad(j, 2)));
      Ad(j, 5) = complex<T>(x, imag(Ad(j, 5)));
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "SwapRealColumn incorrect" << endl;
      abort();
    }

  A.SwapImagColumn(3, 0);
  for (int j = 0; j < m; j++)
    {
      x = imag(Ad(j, 3));
      Ad(j, 3) = complex<T>(real(Ad(j, 3)), imag(Ad(j, 0)));
      Ad(j, 0) = complex<T>(real(Ad(j, 0)), x);
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "SwapImagColumn incorrect" << endl;
      abort();
    }
  
  IVect new_ind(3);
  x = A.ValReal(0, 5);
  new_ind(0) = 1; new_ind(1) = 2; new_ind(2) = 3;
  A.ReplaceRealIndexColumn(5, new_ind);
  Ad(0, 5) = complex<T>(zero, imag(Ad(0, 5)));
  Ad(1, 5) = complex<T>(x, imag(Ad(1, 5)));
  
  new_ind.Reallocate(2);
  new_ind(0) = 1; new_ind(1) = 3;
  x = A.ValImag(0, 1);
  y = A.ValImag(1, 1);
  A.ReplaceImagIndexColumn(1, new_ind);
  Ad(0, 1) = complex<T>(real(Ad(0, 1)), zero);
  Ad(1, 1) = complex<T>(real(Ad(1, 1)), x);
  Ad(3, 1) = complex<T>(real(Ad(3, 1)), y);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReplaceIndexRealRow incorrect" << endl;
      abort();
    }

}

template<class T, class Prop, class Storage, class T2, class Prop2, class Storage2>
void CheckAddMatrix(Matrix<complex<T>, Prop, Storage>& B, Matrix<T2, Prop2, Storage2>& Ad)
{
  T zero;
  SetComplexZero(zero);
  Vector<T> vec(30);
  FillRand(vec);
  
  B.Clear();
  B.Reallocate(4, 7);
  B.AddInteraction(0, 2, complex<T>(vec(0), zero));
  B.AddInteraction(0, 5, complex<T>(vec(1), vec(13)));
  B.AddInteraction(0, 3, complex<T>(zero, vec(2)));
  B.AddInteraction(1, 4, complex<T>(vec(3), vec(14)));
  IVect num(3); Vector<complex<T> > val(3);
  num(0) = 0; num(1) = 2; num(2) = 4;
  val(0) = complex<T>(zero, vec(4)); val(1) = complex<T>(zero, vec(5));
  val(2) = complex<T>(zero, vec(6));
  B.AddInteractionRow(1, 3, num, val);
  num(0) = 1; num(1) = 5;
  val(0) = complex<T>(vec(7), vec(15)); val(1) = complex<T>(vec(8), zero);
  B.AddInteractionRow(2, 2, num, val);
  B.AddInteraction(3, 4, complex<T>(vec(9), vec(16)));
  num(0) = 3; num(1) = 0; num(2) = 2;
  val(0) = complex<T>(vec(10), zero); val(1) = complex<T>(vec(11), zero); val(2) = complex<T>(vec(12), vec(17));
  B.AddInteractionColumn(3, 3, num, val);
  
  Ad.Fill(complex<T>(zero, zero));
  Ad(0, 2) = complex<T>(vec(0), zero); Ad(0, 3) = complex<T>(vec(11), vec(2)); Ad(0, 5) = complex<T>(vec(1), vec(13));
  Ad(1, 4) = complex<T>(vec(3), vec(14)+vec(6)); Ad(1, 0) = complex<T>(zero, vec(4)); Ad(1, 2) = complex<T>(zero, vec(5));
  Ad(2, 1) = complex<T>(vec(7), vec(15)); Ad(2, 3) = complex<T>(vec(12), vec(17)); Ad(2, 5) = complex<T>(vec(8), zero);
  Ad(3, 4) = complex<T>(vec(9), vec(16)); Ad(3, 3) = complex<T>(vec(10), zero);
  
  if (!EqualMatrix(Ad, B))
    {
      cout << "AddInteraction incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckAddMatrix(Matrix<complex<T>, Prop, RowComplexSparse>& B, Matrix<T2, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckAddMatrix(Matrix<complex<T>, Prop, ColComplexSparse>& B, Matrix<T2, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Storage>
void CheckMatrix(Matrix<complex<T>, Prop, Storage>& A)
{
  if (!Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }

  T zero, x, y;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  
  A.Reallocate(4, 7);
  if ( (A.GetM() != 4) || (A.GetN() != 7) || (A.GetDataSize() != 0) || (A.GetRealDataSize() != 0) || (A.GetImagDataSize() != 0) )
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }

  Matrix<complex<T>, Prop, Storage> B(4, 7);
  if ( (B.GetM() != 4) || (B.GetN() != 7) || (B.GetDataSize() != 0) || (B.GetRealDataSize() != 0) || (B.GetImagDataSize() != 0))
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }
  
  Vector<T> vec(30);
  FillRand(vec);
  A.GetReal(0, 0) = vec(0); A.GetReal(0, 2) = vec(1); A.GetReal(0, 5) = vec(2);
  A.GetReal(2, 0) = vec(7); A.GetReal(2, 2) = vec(8); A.GetReal(2, 4) = vec(9); A.GetReal(2, 6) = vec(10);
  A.GetReal(1, 0) = vec(3); A.GetReal(1, 3) = vec(4); A.GetReal(1, 4) = vec(5); A.GetReal(1, 6) = vec(6);
  A.GetReal(3, 5) = vec(11); A.GetReal(3, 1) = vec(12); A.GetReal(3, 4) = vec(13);
  
  A.GetImag(0, 1) = vec(14); A.GetImag(0, 4) = vec(15); A.GetImag(0, 5) = vec(16);
  A.GetImag(3, 2) = vec(17); A.GetImag(3, 5) = vec(18);
  A.GetImag(2, 0) = vec(20); A.GetImag(2, 2) = vec(21); A.GetImag(2, 3) = vec(22);
  A.GetImag(1, 1) = vec(23); A.GetImag(1, 6) = vec(24); A.GetImag(1, 5) = vec(25); A.GetImag(1, 3) = vec(26);

  Matrix<complex<T>, General, RowMajor> Ad(4, 7);
  Ad.Fill(complex<T>(zero, zero));
  Ad(0, 0) = complex<T>(vec(0), zero); Ad(0, 1) = complex<T>(zero, vec(14)); Ad(0, 2) = complex<T>(vec(1), zero);
  Ad(0, 4) = complex<T>(zero, vec(15)); Ad(0, 5) = complex<T>(vec(2), vec(16));
  
  Ad(1, 0) = complex<T>(vec(3), zero); Ad(1, 1) = complex<T>(zero, vec(23)); Ad(1, 3) = complex<T>(vec(4), vec(26));
  Ad(1, 4) = complex<T>(vec(5), zero); Ad(1, 5) = complex<T>(zero, vec(25)); Ad(1, 6) = complex<T>(vec(6), vec(24));
  
  Ad(2, 0) = complex<T>(vec(7), vec(20)); Ad(2, 2) = complex<T>(vec(8), vec(21));
  Ad(2, 3) = complex<T>(zero, vec(22)); Ad(2, 4) = complex<T>(vec(9), zero); Ad(2, 6) = complex<T>(vec(10), zero);
  
  Ad(3, 2) = complex<T>(zero, vec(17)); Ad(3, 5) = complex<T>(vec(11), vec(18));
  Ad(3, 1) = complex<T>(vec(12), zero); Ad(3, 4) = complex<T>(vec(13), zero);
  
  if ( (A.GetRealDataSize() != 14) || (A.GetImagDataSize() != 12) ||  (A.GetDataSize() != 26) )
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "Get incorrect" << endl;
      abort();
    }
  
  A.ValReal(1, 3) = vec(14);
  A.ValReal(2, 6) += vec(15);
  Ad(1, 3) = complex<T>(vec(14), imag(Ad(1, 3)));
  Ad(2, 6) += vec(15);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ValReal incorrect" << endl;
      abort();
    }

  A.ValImag(0, 5) = vec(28);
  A.ValImag(2, 3) += vec(29);
  Ad(0, 5) = complex<T>(real(Ad(0, 5)), vec(28));
  Ad(2, 3) += complex<T>(zero, vec(29));
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ValImag incorrect" << endl;
      abort();
    }
  
  A.Set(3, 6, complex<T>(vec(16), vec(19)));
  Ad.Set(3, 6, complex<T>(vec(16), vec(19)));
  if (!EqualMatrix(A, Ad))
    {
      cout << "Set incorrect" << endl;
      abort();
    }

  Matrix<complex<T>, Prop, Storage> C(A);

  if (!EqualMatrix(A, C))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  Matrix<complex<T>, General, RowMajor> Cd(Ad);  
  
  for (int i = 0; i < C.GetM(); i++)
    for (int j = 0; j < C.GetN(); j++)
      C.Set(i, j, complex<T>(zero, zero));
  
  Cd.Fill(complex<T>(zero, zero));
  if (!EqualMatrix(C, Cd))
    {
      cout << "Set incorrect" << endl;
      abort();
    }

  C.Clear();
  if ( (C.GetM() != 0) || (C.GetN() != 0) || (C.GetDataSize() != 0) || (C.GetRealDataSize() != 0) || (C.GetImagDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
 
  CheckRowMatrix(A, Ad);
  
  Matrix<complex<T>, General, RowMajor> Bd(Ad);  
  C = A;
  if (!EqualMatrix(A, C))
    {
      cout << "Operator = incorrect" << endl;
      abort();
    }
  
  CheckAddMatrix(B, Ad);

  A.WriteText("toto.dat");
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      DISP(A); DISP(B);
      cout << "ReadText/WriteText incorrect" << endl;
      abort();
    }

#ifndef MONTJOIE_WITH_MPFR
  A.Write("totob.dat");
  B.Clear();
  B.Read("totob.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "Read/Write incorrect" << endl;
      abort();
    }
#endif
  
  A.Fill(complex<T>(x, y));
  for (int i = 0; i < Ad.GetM(); i++)
    for (int j = 0; j < Ad.GetN(); j++)
      {
        if (real(Bd(i, j)) != zero)
          Bd(i, j) = complex<T>(x, imag(Bd(i, j)));
        
        if (imag(Bd(i, j)) != zero)
          Bd(i, j) = complex<T>(real(Bd(i, j)), y);
      }
  
  if (!EqualMatrix(A, Bd))
    {
      cout << "Fill incorrect" << endl;
      abort();
    }
  
  A.Reallocate(10, 10);
  Bd.Reallocate(10, 10);
  
  A.SetIdentity();
  Bd.SetIdentity();
  if (!EqualMatrix(A, Bd))
    {
      cout << "SetIdentity incorrect" << endl;
      abort();
    }  
}


template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<complex<T>, Prop, ArrayRowSymComplexSparse>& A,
                    Matrix<complex<T>, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w, x2, y2, z2;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  GetRand(x2); GetRand(y2);
  GetRand(z); GetRand(w);   GetRand(z2);

  int n = Ad.GetN();
  for (int i = 0; i < n; i++)
    {
      int nb_r = 0, nb_i = 0;
      for (int j = i; j < n; j++)
        {
          if (real(Ad(i, j)) != zero)
            nb_r++;
          
          if (imag(Ad(i, j)) != zero)
            nb_i++;
        }
      
      if (nb_r != A.GetRealRowSize(i))
        {
          cout << "GetRealRowSize incorrect" << endl;
          abort();
        }

      if (nb_i != A.GetImagRowSize(i))
        {
          cout << "GetImagRowSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb_r; j++)
        {
          int jcol = A.IndexReal(i, j);
          if ((real(Ad(i, jcol)) != A.ValueReal(i, j)) || isnan(real(Ad(i, jcol))) || isnan(A.ValueReal(i, j)))
            {
              cout << "IndexReal/ValueReal incorrect" << endl;
              abort();
            }
        }

      for (int j = 0; j < nb_i; j++)
        {
          int jcol = A.IndexImag(i, j);
          if ((imag(Ad(i, jcol)) != A.ValueImag(i, j)) || isnan(imag(Ad(i, jcol))) || isnan(A.ValueImag(i, j)))
            {
              cout << "IndexImag/ValueImag incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateRealRow(2, 2);
  A.IndexReal(2, 0) = 2; A.IndexReal(2, 1) = 5; 
  A.ValueReal(2, 0) = x; A.ValueReal(2, 1) = y; 
  for (int j = 2; j < n; j++)
    Ad(2, j) = complex<T>(zero, imag(Ad(2, j)));
  
  Ad(2, 2) = complex<T>(x, imag(Ad(2, 2)));
  Ad(2, 5) = complex<T>(y, imag(Ad(2, 5)));
  
  A.ClearRealRow(4);
  for (int i = 4; i < n; i++)
    Ad(4, i) = complex<T>(zero, imag(Ad(4, i)));

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateRealRow/ClearRealRow incorrect" << endl;
      abort();
    }

  A.ReallocateImagRow(4, 2);
  A.IndexImag(4, 0) = 4; A.IndexImag(4, 1) = 6; 
  A.ValueImag(4, 0) = x2; A.ValueImag(4, 1) = y2; 
  for (int j = 4; j < n; j++)
    Ad(4, j) = complex<T>(real(Ad(4, j)), zero);
  
  Ad(4, 4) = complex<T>(real(Ad(4, 4)), x2);
  Ad(4, 6) = complex<T>(real(Ad(4, 6)), y2);
  
  A.ClearImagRow(5);
  for (int i = 5; i < n; i++)
    Ad(5, i) = complex<T>(real(Ad(5, i)), zero);

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateImagRow/ClearImagRow incorrect" << endl;
      abort();
    }
  
  A.ResizeRealRow(2, 4);
  A.IndexReal(2, 2) = 3; A.IndexReal(2, 3) = 6;
  A.ValueReal(2, 2) = z; A.ValueReal(2, 3) = w;
  A.AssembleRealRow(2);
  Ad(2, 3) = complex<T>(z, imag(Ad(2, 3)));
  Ad(2, 6) = complex<T>(w, imag(Ad(2, 6)));
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeRealRow incorrect" << endl;
      abort();
    }
  
  A.ResizeImagRow(4, 3);
  A.IndexImag(4, 2) = 5;
  A.ValueImag(4, 2) = z2;
  A.AssembleImagRow(4);
  Ad(4, 5) = complex<T>(real(Ad(4, 5)), z2);
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeImagRow incorrect" << endl;
      abort();
    }
}


template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<complex<T>, Prop, RowSymComplexSparse>& A,
                    Matrix<complex<T>, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<complex<T>, Prop, ColSymComplexSparse>& A,
                    Matrix<complex<T>, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<complex<T>, Prop, ArrayColSymComplexSparse>& A,
                    Matrix<complex<T>, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w, x2, y2, z2;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  GetRand(x2); GetRand(y2);
  GetRand(z); GetRand(w); GetRand(z2);

  int m = Ad.GetM();
  for (int i = 0; i < m; i++)
    {
      int nb_r = 0, nb_i = 0;
      for (int j = 0; j <= i; j++)
        {
          if (real(Ad(j, i)) != zero)
            nb_r++;

          if (imag(Ad(j, i)) != zero)
            nb_i++;
        }
      
      if (nb_r != A.GetRealColumnSize(i))
        {
          cout << "GetRealColumnSize incorrect" << endl;
          abort();
        }

      if (nb_i != A.GetImagColumnSize(i))
        {
          cout << "GetImagColumnSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb_r; j++)
        {
          int jrow = A.IndexReal(i, j);
          if ((real(Ad(jrow, i)) != A.ValueReal(i, j)) || isnan(real(Ad(jrow, i))) || isnan(A.ValueReal(i, j)))
            {
              cout << "IndexReal/ValueReal incorrect" << endl;
              abort();
            }
        }

      for (int j = 0; j < nb_i; j++)
        {
          int jrow = A.IndexImag(i, j);
          if ((imag(Ad(jrow, i)) != A.ValueImag(i, j)) || isnan(imag(Ad(jrow, i))) || isnan(A.ValueImag(i, j)))
            {
              cout << "IndexImag/ValueImag incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateRealColumn(2, 2);
  A.IndexReal(2, 0) = 0; A.IndexReal(2, 1) = 1; 
  A.ValueReal(2, 0) = x; A.ValueReal(2, 1) = y; 
  Ad(0, 2) = complex<T>(x, imag(Ad(0, 2)));
  Ad(1, 2) = complex<T>(y, imag(Ad(1, 2)));
  Ad(2, 2) = complex<T>(zero, imag(Ad(2, 2)));
  
  A.ClearRealColumn(1);
  for (int i = 0; i <= 1; i++)
    Ad(i, 1) = complex<T>(zero, imag(Ad(i, 1)));
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateRealColumn/ClearRealColumn incorrect" << endl;
      abort();
    }

  A.ReallocateImagColumn(3, 2);
  A.IndexImag(3, 0) = 0; A.IndexImag(3, 1) = 2; 
  A.ValueImag(3, 0) = x2; A.ValueImag(3, 1) = y2; 
  Ad(0, 3) = complex<T>(real(Ad(0, 3)), x2);
  Ad(1, 3) = complex<T>(real(Ad(1, 3)), zero);
  Ad(2, 3) = complex<T>(real(Ad(2, 3)), y2);
  Ad(3, 3) = complex<T>(real(Ad(3, 3)), zero);
  
  A.ClearImagColumn(2);
  for (int i = 0; i <= 2; i++)
    Ad(i, 2) = complex<T>(real(Ad(i, 2)), zero);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateImagColumn/ClearImagColumn incorrect" << endl;
      abort();
    }
  
  A.ResizeRealColumn(2, 3);
  A.IndexReal(2, 2) = 2;
  A.ValueReal(2, 2) = z;
  A.AssembleRealColumn(2);
  Ad(2, 2) = complex<T>(z, imag(Ad(2, 2)));
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeRealColumn incorrect" << endl;
      abort();
    }  

  A.ResizeImagColumn(3, 4);
  A.IndexImag(3, 2) = 1;   A.IndexImag(3, 3) = 3;
  A.ValueImag(3, 2) = w; A.ValueImag(3, 3) = z2;
  A.AssembleImagColumn(3);
  Ad(1, 3) = complex<T>(real(Ad(1, 3)), w);
  Ad(3, 3) = complex<T>(real(Ad(3, 3)), z2);
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeImagColumn incorrect" << endl;
      abort();
    }  
}

template<class T, class Prop, class Storage, class T2, class Prop2, class Storage2>
void CheckAddMatrixSym(Matrix<complex<T>, Prop, Storage>& B, Matrix<T2, Prop2, Storage2>& Ad)
{
  T zero;
  SetComplexZero(zero);
  Vector<T> vec(30);
  FillRand(vec);
  
  B.Clear();
  B.Reallocate(6, 6);
  B.AddInteraction(0, 2, complex<T>(vec(0), zero));
  B.AddInteraction(0, 5, complex<T>(vec(1), vec(15)));
  B.AddInteraction(0, 3, complex<T>(zero, vec(2)));
  B.AddInteraction(1, 4, complex<T>(vec(3), zero));
  IVect num(3); Vector<complex<T> > val(3);
  num(0) = 0; num(1) = 2; num(2) = 4;
  val(0) = complex<T>(zero, vec(4)); val(1) = complex<T>(zero, vec(5));
  val(2) = complex<T>(vec(6), vec(16));
  B.AddInteractionRow(1, 3, num, val);
  num(0) = 1; num(1) = 5;
  val(0) = complex<T>(vec(7), zero); val(1) = complex<T>(vec(8), zero);
  B.AddInteractionRow(2, 2, num, val);
  B.AddInteraction(3, 4, complex<T>(vec(9), vec(17)));
  num(0) = 5; num(1) = 0; num(2) = 2;
  val(0) = complex<T>(zero, vec(10));
  val(1) = complex<T>(vec(11), zero); val(2) = complex<T>(vec(12), vec(18));
  B.AddInteractionColumn(3, 3, num, val);
  B.AddInteraction(4, 4, complex<T>(vec(13), vec(19)));
  
  Ad.Reallocate(6, 6);
  Ad.Fill(0);
  Ad(0, 2) = complex<T>(vec(0), zero); Ad(0, 3) = complex<T>(vec(11), vec(2));
  Ad(0, 5) = complex<T>(vec(1), vec(15));;
  Ad(1, 4) = complex<T>(vec(3)+vec(6), vec(16)); Ad(1, 2) = complex<T>(zero, vec(5));
  Ad(2, 3) = complex<T>(vec(12), vec(18)); Ad(2, 5) = complex<T>(vec(8), zero);
  Ad(3, 4) = complex<T>(vec(9), vec(17)); Ad(4, 4) = complex<T>(vec(13), vec(19));
  
  if (!EqualMatrix(Ad, B))
    {
      cout << "AddInteraction incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckAddMatrixSym(Matrix<complex<T>, Prop, RowSymComplexSparse>& B,
                       Matrix<T2, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class T2, class Prop2, class Storage2>
void CheckAddMatrixSym(Matrix<complex<T>, Prop, ColSymComplexSparse>& B,
                       Matrix<T2, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Storage>
void CheckMatrixSym(Matrix<complex<T>, Prop, Storage>& A)
{
  if (!Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }

  T zero, x, x2;
  SetComplexZero(zero);
  GetRand(x); GetRand(x2);
  
  A.Reallocate(7, 7);
  if ( (A.GetM() != 7) || (A.GetN() != 7) || (A.GetDataSize() != 0))
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }
  
  Matrix<complex<T>, Prop, Storage> B(7, 7);
  if ( (B.GetM() != 7) || (B.GetN() != 7) || (B.GetDataSize() != 0))
    {
      cout << "Constructor incorrect " << endl;
      abort();
    }
  
  Vector<T> vec(40);
  FillRand(vec);
  A.GetReal(0, 0) = vec(0); A.GetReal(0, 2) = vec(1); A.GetReal(0, 5) = vec(2);
  A.GetReal(2, 2) = vec(6); A.GetReal(2, 4) = vec(7); A.GetReal(2, 6) = vec(16);
  A.GetReal(1, 1) = vec(3); A.GetReal(1, 3) = vec(4); A.GetReal(1, 4) = vec(5);
  A.GetReal(3, 5) = vec(8); A.GetReal(3, 3) = vec(9); A.GetReal(3, 4) = vec(10);
  A.GetReal(4, 4) = vec(11);  A.GetReal(4, 5) = vec(12);
  A.GetReal(5, 5) = vec(13); A.GetReal(5, 6) = vec(14);
  A.GetReal(6, 6) = vec(15);

  A.GetImag(0, 0) = vec(17); A.GetImag(0, 1) = vec(18); A.GetImag(0, 4) = vec(19); A.GetImag(0, 6) = vec(20);
  A.GetImag(2, 3) = vec(21); A.GetImag(2, 4) = vec(22);
  A.GetImag(1, 1) = vec(23); A.GetImag(1, 3) = vec(24); A.GetImag(1, 5) = vec(25);
  A.GetImag(3, 6) = vec(26); A.GetImag(3, 3) = vec(27); A.GetImag(3, 4) = vec(28);
  A.GetImag(4, 4) = vec(29);  A.GetImag(4, 6) = vec(30);
  A.GetImag(5, 5) = vec(31); A.GetImag(6, 6) = vec(32);
  
  Matrix<complex<T>, Symmetric, RowSymPacked> Ad(7, 7);
  Ad.Fill(complex<T>(zero, zero));
  Ad(0, 0) = complex<T>(vec(0), vec(17)); Ad(0, 1) = complex<T>(zero, vec(18));
  Ad(0, 2) = complex<T>(vec(1), zero); Ad(0, 4) = complex<T>(zero, vec(19));
  Ad(0, 5) = complex<T>(vec(2), zero); Ad(0, 6) = complex<T>(zero, vec(20));
  
  Ad(2, 2) = complex<T>(vec(6), zero); Ad(2, 3) = complex<T>(zero, vec(21));
  Ad(2, 4) = complex<T>(vec(7), vec(22)); Ad(2, 6) = complex<T>(vec(16), zero);
  
  Ad(1, 1) = complex<T>(vec(3), vec(23)); Ad(1, 3) = complex<T>(vec(4), vec(24));
  Ad(1, 4) = complex<T>(vec(5), zero); Ad(1, 5) = complex<T>(zero, vec(25));
  
  Ad(3, 3) = complex<T>(vec(9), vec(27)); Ad(3, 4) = complex<T>(vec(10), vec(28));
  Ad(3, 5) = complex<T>(vec(8), zero);  Ad(3, 6) = complex<T>(zero, vec(26)); 
  
  Ad(4, 4) = complex<T>(vec(11), vec(29));  Ad(4, 5) = complex<T>(vec(12), zero); Ad(4, 6) = complex<T>(zero, vec(30));
  Ad(5, 5) = complex<T>(vec(13), vec(31)); Ad(5, 6) = complex<T>(vec(14), zero);
  Ad(6, 6) = complex<T>(vec(15), vec(32));
  
  //DISP(A);
  if ( (A.GetImagDataSize() != 16) || (A.GetRealDataSize() != 17) || (A.GetDataSize() != 33))
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "GetReal/GetImag incorrect" << endl;
      abort();
    }
  
  A.ValReal(1, 3) = vec(33);
  A.ValReal(2, 4) += vec(34);
  Ad(1, 3) = complex<T>(vec(33), imag(Ad(1, 3)));
  Ad(2, 4) += vec(34);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ValReal incorrect" << endl;
      abort();
    }

  A.ValImag(1, 1) = vec(35);
  A.ValImag(2, 4) += vec(36);
  Ad(1, 1) = complex<T>(real(Ad(1, 1)), vec(35));
  Ad(2, 4) += complex<T>(zero, vec(36));
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ValImag incorrect" << endl;
      abort();
    }
  
  A.Set(3, 5, complex<T>(vec(39), vec(38)));
  Ad.Set(3, 5, complex<T>(vec(39), vec(38)));
  if (!EqualMatrix(A, Ad))
    {
      cout << "Set incorrect" << endl;
      abort();
    }

  Matrix<complex<T>, Prop, Storage> C(A);

  if (!EqualMatrix(A, C))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  Matrix<complex<T>, Symmetric, RowSymPacked> Cd(Ad);  
  
  for (int i = 0; i < C.GetM(); i++)
    for (int j = 0; j < C.GetN(); j++)
      C.Set(i, j, complex<T>(zero, zero));
  
  Cd.Fill(complex<T>(zero, zero));
  if (!EqualMatrix(C, Cd))
    {
      cout << "Set incorrect" << endl;
      abort();
    }

  C.Clear();
  if ( (C.GetM() != 0) || (C.GetN() != 0) || (C.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
 
  CheckRowMatrix(A, Ad);
  
  Matrix<complex<T>, Symmetric, RowSymPacked> Bd(Ad);  
  C = A;
  if (!EqualMatrix(A, C))
    {
      cout << "Operator = incorrect" << endl;
      abort();
    }
  
  CheckAddMatrixSym(B, Ad);
  
  A.WriteText("toto.dat");
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "ReadText/WriteText incorrect" << endl;
      abort();
    }

#ifndef MONTJOIE_WITH_MPFR
  A.Write("totob.dat");
  B.Clear();
  B.Read("totob.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "Read/Write incorrect" << endl;
      abort();
    }
#endif
  
  A.Fill(complex<T>(x, x2));
  for (int i = 0; i < Bd.GetM(); i++)
    for (int j = 0; j < Bd.GetN(); j++)
      {
      if (real(Bd(i, j)) != zero)
          Bd(i, j) = complex<T>(x, imag(Bd(i, j)));
        
        if (imag(Bd(i, j)) != zero)
          Bd(i, j) = complex<T>(real(Bd(i, j)), x2);
      }
  
  if (!EqualMatrix(A, Bd))
    {
      cout << "Fill incorrect" << endl;
      abort();
    }
  
  A.Reallocate(10, 10);
  Bd.Reallocate(10, 10);
  
  A.SetIdentity();
  Bd.SetIdentity();
  if (!EqualMatrix(A, Bd))
    {
      cout << "SetIdentity incorrect" << endl;
      abort();
    }  
}

int main(int argc, char** argv)
{
  threshold = 1e-12;
  cout.precision(15);
  
  {
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    CheckMatrix(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
    CheckMatrix(A);
  }
#endif

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    CheckMatrixSym(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
    CheckMatrixSym(A);
  }
#endif

  {
    Matrix<Complex_wp, General, RowComplexSparse> A;
    CheckMatrix(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, General, ColComplexSparse> A;
    CheckMatrix(A);
  }
#endif

  {
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    CheckMatrixSym(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
    CheckMatrixSym(A);
  }
#endif

  cout << "All tests passed successfully" << endl;
  
  return 0;
}
