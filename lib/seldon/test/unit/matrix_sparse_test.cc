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

Real_wp threshold = 1e-12;

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void FillRand(Vector<T> & x)
{
  typedef typename ClassComplexType<T>::Treal Treal;
  x.FillRand();
  Mlt(Treal(1e-9), x);
}

template<class T>
void FillRand(Vector<complex<T> > & x)
{
  for (int i = 0; i < x.GetM(); i++)
    x(i) = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T, class Prop1, class Storage1, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop1, Storage1>& A,
		 const Matrix<T, Prop2, Storage2>& B)
{
  if ( (A.GetM() != B.GetM()) || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if (( abs(A(i, j) - B(i, j) ) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          //abort();
          return false;
        }
  
  return true;
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, RowSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)
{}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ColSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)
{}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ArrayRowSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  GetRand(z); GetRand(w);

  int n = Ad.GetN();
  for (int i = 0; i < Ad.GetM(); i++)
    {
      int nb = 0;
      for (int j = 0; j < n; j++)
        if (Ad(i, j) != zero)
          nb++;
      
      if (nb != A.GetRowSize(i))
        {
          cout << "GetRowSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb; j++)
        {
          int jcol = A.Index(i, j);
          if ((Ad(i, jcol) != A.Value(i, j)) || isnan(Ad(i, jcol)) || isnan(A.Value(i, j)))
            {
              cout << "Index/Value incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateRow(2, 2);
  A.Index(2, 0) = 2; A.Index(2, 1) = 5; 
  A.Value(2, 0) = x; A.Value(2, 1) = y; 
  Ad(2,0) = zero; Ad(2,1) = zero; Ad(2,2) = x;
  Ad(2,3) = zero; Ad(2,4) = zero; Ad(2,5) = y; Ad(2,6) = zero;
  
  A.ClearRow(3);
  for (int i = 0; i < n; i++)
    Ad(3, i) = zero;

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateRow/ClearRow incorrect" << endl;
      abort();
    }
  
  A.ResizeRow(2, 4);
  A.Index(2, 2) = 1; A.Index(2, 3) = 4;
  A.Value(2, 2) = z; A.Value(2, 3) = w;
  A.AssembleRow(2);
  Ad(2, 1) = z; Ad(2, 4) = w;
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeRow incorrect" << endl;
      abort();
    }
  
  A.SwapRow(1, 3);
  for (int j = 0; j < n; j++)
    {
      x = Ad(1, j);
      Ad(1, j) = Ad(3, j);
      Ad(3, j) = x;
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "SwapRow incorrect" << endl;
      abort();
    }
  
  IVect new_ind(3);
  new_ind(0) = 1; new_ind(1) = 4; new_ind(2) = 5;
  A.ReplaceIndexRow(0, new_ind);
  x = Ad(0, 0); y = Ad(0, 2);
  Ad(0, 0) = zero; Ad(0, 1) = x; Ad(0, 2) = zero; Ad(0, 4) = y;
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReplaceIndexRow incorrect" << endl;
      abort();
    }

}


template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ArrayColSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  GetRand(z); GetRand(w);

  int m = Ad.GetM();
  for (int i = 0; i < Ad.GetN(); i++)
    {
      int nb = 0;
      for (int j = 0; j < m; j++)
        if (Ad(j, i) != zero)
          nb++;
      
      if (nb != A.GetColumnSize(i))
        {
          cout << "GetColumnSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb; j++)
        {
          int jrow = A.Index(i, j);
          if ((Ad(jrow, i) != A.Value(i, j)) || isnan(Ad(jrow, i)) || isnan(A.Value(i, j)))
            {
              cout << "Index/Value incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateColumn(2, 2);
  A.Index(2, 0) = 2; A.Index(2, 1) = 3; 
  A.Value(2, 0) = x; A.Value(2, 1) = y; 
  Ad(0, 2) = zero; Ad(1, 2) = zero; Ad(2, 2) = x;
  Ad(3, 2) = y;
  
  A.ClearColumn(5);
  for (int i = 0; i < m; i++)
    Ad(i, 5) = zero;

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateColumn/ClearColumn incorrect" << endl;
      abort();
    }
  
  A.ResizeColumn(1, 3);
  A.Index(1, 1) = 0; A.Index(1, 2) = 2;
  A.Value(1, 1) = z; A.Value(1, 2) = w;
  A.AssembleColumn(1);
  Ad(0, 1) = z; Ad(2, 1) = w;
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeColumn incorrect" << endl;
      abort();
    }
  
  A.SwapColumn(2, 5);
  for (int j = 0; j < m; j++)
    {
      x = Ad(j, 2);
      Ad(j, 2) = Ad(j, 5);
      Ad(j, 5) = x;
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "SwapColumn incorrect" << endl;
      abort();
    }
  
  IVect new_ind(2);
  new_ind(0) = 0; new_ind(1) = 2;
  A.ReplaceIndexColumn(5, new_ind);
  x = Ad(2, 5); y = Ad(3, 5);
  Ad(0, 5) = x; Ad(1, 5) = zero; Ad(2, 5) = y; Ad(3, 5) = zero;
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "ReplaceIndexRow incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Storage, class Prop2, class Storage2>
void CheckAddMatrix(Matrix<T, Prop, Storage>& B,
                    Matrix<T, Prop2, Storage2>& Ad)
{
  T zero; SetComplexZero(zero);
  Vector<T> vec(20);
  FillRand(vec);

  B.Clear();
  B.Reallocate(4, 7);
  B.AddInteraction(0, 2, vec(0));
  B.AddInteraction(0, 5, vec(1));
  B.AddInteraction(0, 3, vec(2));
  B.AddInteraction(1, 4, vec(3));
  IVect num(3); Vector<T> val(3);
  num(0) = 0; num(1) = 2; num(2) = 4;
  val(0) = vec(4); val(1) = vec(5); val(2) = vec(6);
  B.AddInteractionRow(1, 3, num, val);
  num(0) = 1; num(1) = 5;
  val(0) = vec(7); val(1) = vec(8);
  B.AddInteractionRow(2, 2, num, val);
  B.AddInteraction(3, 4, vec(9));
  num(0) = 3; num(1) = 0; num(2) = 2;
  val(0) = vec(10); val(1) = vec(11); val(2) = vec(12);
  B.AddInteractionColumn(3, 3, num, val);
  
  Ad.Fill(zero);
  Ad(0, 2) = vec(0); Ad(0, 3) = vec(2)+vec(11); Ad(0, 5) = vec(1);
  Ad(1, 4) = vec(3)+vec(6); Ad(1, 0) = vec(4); Ad(1, 2) = vec(5);
  Ad(2, 1) = vec(7); Ad(2, 3) = vec(12); Ad(2, 5) = vec(8);
  Ad(3, 4) = vec(9); Ad(3, 3) = vec(10);
  
  if (!EqualMatrix(Ad, B))
    {
      cout << "AddInteraction incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckAddMatrix(Matrix<T, Prop, RowSparse>& B,
                    Matrix<T, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckAddMatrix(Matrix<T, Prop, ColSparse>& B,
                    Matrix<T, Prop2, Storage2>& Ad)
{
}
  
template<class T, class Prop, class Storage>
void CheckMatrix(Matrix<T, Prop, Storage>& A)
{
  if (!Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }

  T zero, x;
  SetComplexZero(zero);
  GetRand(x);
  
  A.Reallocate(4, 7);
  if ( (A.GetM() != 4) || (A.GetN() != 7) || (A.GetDataSize() != 0))
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }

  Matrix<T, Prop, Storage> B(4, 7);
  if ( (B.GetM() != 4) || (B.GetN() != 7) || (B.GetDataSize() != 0))
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }
  
  Vector<T> vec(20);
  FillRand(vec);
  A.Get(0, 0) = vec(0); A.Get(0, 2) = vec(1); A.Get(0, 5) = vec(2);
  A.Get(2, 0) = vec(7); A.Get(2, 2) = vec(8); A.Get(2, 4) = vec(9); A.Get(2, 6) = vec(10);
  A.Get(1, 0) = vec(3); A.Get(1, 3) = vec(4); A.Get(1, 4) = vec(5); A.Get(1, 6) = vec(6);
  A.Get(3, 5) = vec(11); A.Get(3, 1) = vec(12); A.Get(3, 4) = vec(13);
  
  Matrix<T, General, RowMajor> Ad(4, 7);
  Ad.Fill(zero);
  Ad.Get(0, 0) = vec(0); Ad.Get(0, 2) = vec(1); Ad.Get(0, 5) = vec(2);
  Ad.Get(1, 0) = vec(3); Ad.Get(1, 3) = vec(4); Ad.Get(1, 4) = vec(5); Ad.Get(1, 6) = vec(6);
  Ad.Get(2, 0) = vec(7); Ad.Get(2, 2) = vec(8); Ad.Get(2, 4) = vec(9); Ad.Get(2, 6) = vec(10);
  Ad.Get(3, 5) = vec(11); Ad.Get(3, 1) = vec(12); Ad.Get(3, 4) = vec(13);
  
  if ( (A.GetDataSize() != 14) || (A.GetNonZeros() != 14))
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "Get incorrect" << endl;
      abort();
    }
  
  A.Val(1, 3) = vec(14);
  A.Val(2, 6) += vec(15);
  Ad.Val(1, 3) = vec(14);
  Ad.Val(2, 6) += vec(15);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "Val incorrect" << endl;
      abort();
    }
  
  A.Set(3, 6, vec(16));
  Ad.Set(3, 6, vec(16));
  if (!EqualMatrix(A, Ad))
    {
      cout << "Set incorrect" << endl;
      abort();
    }

  Matrix<T, Prop, Storage> C(A);

  if (!EqualMatrix(A, C))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  C.Clear();
  if ( (C.GetM() != 0) || (C.GetN() != 0) || (C.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
 
  CheckRowMatrix(A, Ad);
  
  Matrix<T, General, RowMajor> Bd(Ad);  
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
  
  A.Fill(x);
  for (int i = 0; i < Ad.GetM(); i++)
    for (int j = 0; j < Ad.GetN(); j++)
      {
        if (Bd(i, j) != zero)
          Bd(i, j) = x;
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

  A.Reallocate(120, 50);
  Bd.Reallocate(120, 50);
  Bd.Zero();
  for (int n = 0; n < 2000; n++)
    {
      int i = rand()%A.GetM();
      int j = rand()%A.GetN();
      GetRand(x);
      Bd(i, j) += x;
      A.Get(i, j) += x;
    }
  
  if (!EqualMatrix(A, Bd))
    {
      cout << "Get incorrect" << endl;
      abort();
    }  
}


template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ArrayRowSymSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  GetRand(z); GetRand(w);

  int n = Ad.GetN();
  for (int i = 0; i < n; i++)
    {
      int nb = 0;
      for (int j = i; j < n; j++)
        if (Ad(i, j) != zero)
          nb++;
      
      if (nb != A.GetRowSize(i))
        {
          cout << "GetRowSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb; j++)
        {
          int jcol = A.Index(i, j);
          if ((Ad(i, jcol) != A.Value(i, j)) || isnan(Ad(i, jcol)) || isnan(A.Value(i, j)))
            {
              cout << "Index/Value incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateRow(2, 2);
  A.Index(2, 0) = 2; A.Index(2, 1) = 5; 
  A.Value(2, 0) = x; A.Value(2, 1) = y; 
  Ad(2,2) = x; Ad(2,3) = zero; Ad(2,4) = zero; Ad(2,5) = y;
  
  A.ClearRow(3);
  for (int i = 3; i < n; i++)
    Ad(3, i) = zero;

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateRow/ClearRow incorrect" << endl;
      abort();
    }
  
  A.ResizeRow(2, 4);
  A.Index(2, 2) = 3; A.Index(2, 3) = 4;
  A.Value(2, 2) = z; A.Value(2, 3) = w;
  A.AssembleRow(2);
  Ad(2, 3) = z; Ad(2, 4) = w;
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeRow incorrect" << endl;
      abort();
    }
  
}


template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ArrayColSymSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)

{
  T zero, x, y, z, w;
  SetComplexZero(zero);
  GetRand(x); GetRand(y);
  GetRand(z); GetRand(w);

  int m = Ad.GetM();
  for (int i = 0; i < m; i++)
    {
      int nb = 0;
      for (int j = 0; j <= i; j++)
        if (Ad(j, i) != zero)
          nb++;
      
      if (nb != A.GetColumnSize(i))
        {
          cout << "GetColumnSize incorrect" << endl;
          abort();
        }
      
      for (int j = 0; j < nb; j++)
        {
          int jrow = A.Index(i, j);
          if ((Ad(jrow, i) != A.Value(i, j)) || isnan(Ad(jrow, i)) || isnan(A.Value(i, j)))
            {
              cout << "Index/Value incorrect" << endl;
              abort();
            }
        }
    }
  
  A.ReallocateColumn(2, 2);
  A.Index(2, 0) = 0; A.Index(2, 1) = 1; 
  A.Value(2, 0) = x; A.Value(2, 1) = y; 
  Ad(0, 2) = x; Ad(1, 2) = y; Ad(2, 2) = zero;
  
  A.ClearColumn(4);
  for (int i = 0; i <= 4; i++)
    Ad(i, 4) = zero;

  if (!EqualMatrix(A, Ad))
    {
      cout << "ReallocateColumn/ClearColumn incorrect" << endl;
      abort();
    }
  
  A.ResizeColumn(5, 6);
  A.Index(5, 4) = 1; A.Index(5, 5) = 2;
  A.Value(5, 4) = z; A.Value(5, 5) = w;
  A.AssembleColumn(5);
  Ad(1, 5) = z; Ad(2, 5) = w;
  if (!EqualMatrix(A, Ad))
    {
      cout << "ResizeColumn incorrect" << endl;
      abort();
    }  
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, RowSymSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)
{}

template<class T, class Prop, class Prop2, class Storage2>
void CheckRowMatrix(Matrix<T, Prop, ColSymSparse>& A,
                    Matrix<T, Prop2, Storage2>& Ad)
{}


template<class T, class Prop, class Storage, class Prop2, class Storage2>
void CheckAddMatrixSym(Matrix<T, Prop, Storage>& B,
                       Matrix<T, Prop2, Storage2>& Ad)
{
  T zero; SetComplexZero(zero);
  Vector<T> vec(20);
  FillRand(vec);

  B.Clear();
  B.Reallocate(6, 6);
  B.AddInteraction(0, 2, vec(0));
  B.AddInteraction(0, 5, vec(1));
  B.AddInteraction(0, 3, vec(2));
  B.AddInteraction(1, 4, vec(3));
  IVect num(3); Vector<T> val(3);
  num(0) = 0; num(1) = 2; num(2) = 4;
  val(0) = vec(4); val(1) = vec(5); val(2) = vec(6);
  B.AddInteractionRow(1, 3, num, val);
  num(0) = 1; num(1) = 5;
  val(0) = vec(7); val(1) = vec(8);
  B.AddInteractionRow(2, 2, num, val);
  B.AddInteraction(3, 4, vec(9));
  num(0) = 5; num(1) = 0; num(2) = 2;
  val(0) = vec(10); val(1) = vec(11); val(2) = vec(12);
  B.AddInteractionColumn(3, 3, num, val);
  B.AddInteraction(4, 4, vec(13));
  
  Ad.Fill(zero);
  Ad(0, 2) = vec(0); Ad(0, 3) = vec(2)+vec(11); Ad(0, 5) = vec(1);
  Ad(1, 4) = vec(3)+vec(6); Ad(1, 2) = vec(5);
  Ad(2, 3) = vec(12); Ad(2, 5) = vec(8);
  Ad(3, 4) = vec(9); Ad(4, 4) = vec(13);
  
  if (!EqualMatrix(Ad, B))
    {
      cout << "AddInteraction incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckAddMatrixSym(Matrix<T, Prop, RowSymSparse>& B,
                       Matrix<T, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Prop2, class Storage2>
void CheckAddMatrixSym(Matrix<T, Prop, ColSymSparse>& B,
                       Matrix<T, Prop2, Storage2>& Ad)
{
}

template<class T, class Prop, class Storage>
void CheckMatrixSym(Matrix<T, Prop, Storage>& A)
{
  if (!Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }

  T zero, x;
  SetComplexZero(zero);
  GetRand(x);
  
  A.Reallocate(6, 6);
  if ( (A.GetM() != 6) || (A.GetN() != 6) || (A.GetDataSize() != 0))
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }
  
  Matrix<T, Prop, Storage> B(6, 6);
  if ( (B.GetM() != 6) || (B.GetN() != 6) || (B.GetDataSize() != 0))
    {
      cout << "Reallocate incorrect " << endl;
      abort();
    }
  
  Vector<T> vec(20);
  FillRand(vec);
  A.Get(0, 0) = vec(0); A.Get(0, 2) = vec(1); A.Get(0, 5) = vec(2);
  A.Get(2, 2) = vec(6); A.Get(2, 4) = vec(7);
  A.Get(1, 1) = vec(3); A.Get(1, 3) = vec(4); A.Get(1, 4) = vec(5);
  A.Get(3, 5) = vec(8); A.Get(3, 3) = vec(9); A.Get(3, 4) = vec(10);
  A.Get(4, 4) = vec(11);  A.Get(4, 5) = vec(12);  A.Get(5, 5) = vec(13);
  
  Matrix<T, Symmetric, RowSymPacked> Ad(6, 6);
  Ad.Fill(zero);
  Ad.Get(0, 0) = vec(0); Ad.Get(0, 2) = vec(1); Ad.Get(0, 5) = vec(2);
  Ad.Get(2, 2) = vec(6); Ad.Get(2, 4) = vec(7);
  Ad.Get(1, 1) = vec(3); Ad.Get(1, 3) = vec(4); Ad.Get(1, 4) = vec(5);
  Ad.Get(3, 5) = vec(8); Ad.Get(3, 3) = vec(9); Ad.Get(3, 4) = vec(10);
  Ad.Get(4, 4) = vec(11);  Ad.Get(4, 5) = vec(12);  Ad.Get(5, 5) = vec(13);
  
  if ( (A.GetDataSize() != 14) || (A.GetNonZeros() != 14))
    {
      cout << "GetDataSize incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "Get incorrect" << endl;
      abort();
    }
  
  A.Val(1, 3) = vec(14);
  A.Val(2, 4) += vec(15);
  Ad.Val(1, 3) = vec(14);
  Ad.Val(2, 4) += vec(15);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "Val incorrect" << endl;
      abort();
    }
  
  A.Set(3, 5, vec(16));
  Ad.Set(3, 5, vec(16));
  if (!EqualMatrix(A, Ad))
    {
      cout << "Set incorrect" << endl;
      abort();
    }

  Matrix<T, Prop, Storage> C(A);

  if (!EqualMatrix(A, C))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  C.Clear();
  if ( (C.GetM() != 0) || (C.GetN() != 0) || (C.GetDataSize() != 0))
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
 
  CheckRowMatrix(A, Ad);
  
  Matrix<T, Symmetric, RowSymPacked> Bd(Ad);  
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
  
  A.Fill(x);
  for (int i = 0; i < Ad.GetM(); i++)
    for (int j = 0; j < Ad.GetN(); j++)
      {
        if (Bd(i, j) != zero)
          Bd(i, j) = x;
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
  //srand(time(NULL));

  cout.precision(15);
  
  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckMatrix(A);
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckMatrix(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, General, ArrayColSparse> A;
    CheckMatrix(A);
  }

  {
    Matrix<Complex_wp, General, ArrayColSparse> A;
    CheckMatrix(A);
  }
#endif

  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckMatrixSym(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckMatrixSym(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    CheckMatrixSym(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayColSymSparse> A;
    CheckMatrixSym(A);
  }
#endif
  
  {
    Matrix<Real_wp, General, RowSparse> A, B;
    CheckMatrix(A);
    
    A.Reallocate(30, 40);
    A.FillRand(200);
    A.WriteText("toto.dat");
    B.ReadText("toto.dat");
    if (!EqualMatrix(A, B))
      {
        cout << "FillRand incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, RowSparse> A;
    CheckMatrix(A);
  }
  
  {
    Matrix<Real_wp, General, ColSparse> A, B;
    CheckMatrix(A);
    
    A.Reallocate(30, 40);
    A.FillRand(200);
    A.WriteText("toto.dat");
    B.ReadText("toto.dat");
    if (!EqualMatrix(A, B))
      {
        cout << "FillRand incorrect" << endl;
        abort();
      }
  }

  {
    Matrix<Complex_wp, General, ColSparse> A;
    CheckMatrix(A);
  }

  {
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    CheckMatrixSym(A);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    CheckMatrixSym(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    CheckMatrixSym(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymSparse> A;
    CheckMatrixSym(A);
  }
#endif
  
  std::remove("toto.dat");
  std::remove("totob.dat");
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
