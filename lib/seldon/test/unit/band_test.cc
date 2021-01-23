#include "SeldonLib.hxx"

using namespace Seldon;

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
void CheckBandMatrix(Matrix<T, General, BandedCol>& A, int n, int kl, int ku)
{
  A.Reallocate(n, n, kl, ku);
  
  // checking GetM/GetN
  if ( (A.GetM() != n) || (A.GetN() != n) )
    {
      cout << "GetM/GetN incorrect" << endl;
      abort();
    }
  
  // checking GetKL, GetKU
  if ( (A.GetKL() != kl) || (A.GetKU() != ku) )
    {
      cout << "GetKL/GetKU incorrect" << endl;
      abort();
    }
  
  // checking GetDataSize()
  if (A.GetDataSize() <= n*(kl+ku+1))
    {
      cout << "GetDataSize() incorrect " << endl;
      abort();
    }
    
  // checking Clear
  A.Clear();
  
  if ( (A.GetM() != 0) || (A.GetDataSize() != 0) )
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  // checking Zero
  A.Reallocate(n, n, kl, ku);

  T zero; SetComplexZero(zero);
  T one; SetComplexOne(one);
  typedef typename ClassComplexType<T>::Treal Treal;
  
  A.Zero();
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j)) > threshold) || isnan(A(i, j)))
        {
          DISP(i); DISP(j); DISP(A(i, j));
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  // checking Get/Val
  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand();
  Mlt(Treal(1e-9), B);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if ((i <= j+kl) && (i >=j-ku))
	  {
	    A.Get(i, j) = B(i, j);
	  }
	else
	  B(i, j) = zero;
      }
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(i, j)) || isnan(B(i, j)))
        {
          cout << "Get incorrect" << endl;
          abort();
        }
  
  A.Fill(zero);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((i <= j+kl) && (i >=j-ku))
	A.Val(i, j) = B(i, j);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(i, j)) || isnan(A(i, j)))
        {
          cout << "Val incorrect" << endl;
          abort();
        }
  
  // checking AddInteraction
  Matrix<T, General, RowMajor> C(n, n);
  C.FillRand(); Mlt(Treal(1e-9), C);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if ((i <= j+kl) && (i >=j-ku))
	  A.AddInteraction(i, j, C(i, j));
	else
	  C(i, j) = zero;
      }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j) - C(i, j)) > threshold) || isnan(C(i, j)) || isnan(A(i, j)))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j)); DISP(C(i, j));
          DISP(A(i, j) - B(i, j) - C(i, j));
          cout << "AddInteraction incorrect" << endl;
          abort();
        }

  // checking AddInteractionRow
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((i <= j+kl) && (i >=j-ku))
	A.Val(i, j) = B(i, j);
  
  IVect num(kl+ku+1); Vector<T> val(kl+ku+1);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	int nb_val = 0;
	if ((i <= j+kl) && (i >=j-ku))
	  {
	    num(nb_val) = j;
	    val(nb_val) = C(i, j);
	    nb_val++;
	  }
	
	A.AddInteractionRow(i, nb_val, num, val);
      }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j) - C(i, j)) > threshold) || isnan(C(i, j)) || isnan(A(i, j)))
        {
          cout << "AddInteractionRow incorrect" << endl;
          abort();
        }
  
  // checking ClearRow and operator *=
  A.ClearRow(5);
  T coef_mlt;
  SetComplexReal(to_num<Treal>("1.3"), coef_mlt);
  A *= coef_mlt;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	T val_ref = coef_mlt*(B(i, j) + C(i, j));
	if (i == 5)
	  val_ref = zero;
	
	if ((abs(A(i, j) - val_ref) > threshold) || isnan(A(i, j)) || isnan(val_ref))
	  {
            DISP(i); DISP(val_ref); DISP(A(i, j));
	    cout << "ClearRow incorrect" << endl;
	    abort();
	  }
      }
  
  // checking SetIdentity
  A.SetIdentity();
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        T val_ex; SetComplexZero(val_ex);
        if (i == j)
          SetComplexOne(val_ex);
        
        if ((abs(A(i, j)-val_ex) > threshold) || isnan(A(i, j)))
          {	    
            cout << "SetIdentity incorrect" << endl;
            abort();
          }
      }
  
  // checking Fill
  SetComplexReal(to_num<Treal>("2.5"), coef_mlt);
  A.Fill(coef_mlt);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((i <= j+kl) && (i >=j-ku))
	if ((abs(A(i, j) - coef_mlt) > threshold) || isnan(A(i, j)))
	  {
	    cout << "Fill incorrect" << endl;
	    abort();
	  }
  
  // checking FillRand, Write, WriteText
  A.FillRand(); A *= Treal(1e-9);
  
  A.Write("band_bin.dat");
  
  A.WriteText("band.dat");
  
  // checking conversion from sparse matrices to band matrices
  Matrix<T, General, RowSparse> Acsr(n, n);
  Matrix<T, General, ArrayRowSparse> Asp;
  Acsr.FillRand(80); Mlt(Treal(1e-9), Acsr);
  Acsr.WriteText("mat_csr.dat");  
  Copy(Acsr, Asp);
  Copy(Asp, A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((Acsr(i, j) != A(i, j)) || isnan(Acsr(i, j)) || isnan(A(i, j)))
	{
          DISP(i); DISP(j); DISP(Acsr(i, j)); DISP(A(i, j));
	  cout << "Copy incorrect" << endl;
	  abort();
	}
  
  // checking Add
  Matrix<T, General, BandedCol> Cband;
  Cband.Reallocate(n, n, kl, ku);
  A.Reallocate(n, n, kl, ku);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((i <= j+kl) && (i >=j-ku))
	{
	  A.Get(i, j) = B(i, j);
	  Cband.Get(i, j) = C(i, j);
	}
  
  SetComplexReal(to_num<Treal>("0.9"), coef_mlt);
  Add(coef_mlt, Cband, A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j) - coef_mlt*C(i, j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
  
  // checking Mlt
  A.Reallocate(n, n, kl, ku);
  A.Fill(0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((i <= j+kl) && (i >=j-ku))
	A.Get(i, j) = B(i, j);

  SetComplexReal(to_num<Treal>("3.2"), coef_mlt);
  Mlt(coef_mlt, A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - coef_mlt*B(i, j)) > threshold) || isnan(A(i, j)))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  

  Mlt(one/coef_mlt, A);
  Vector<T> x(n), y(n), y2(n);
  x.FillRand(); x *= 1e-9;
  Mlt(A, x, y);
  Mlt(B, x, y2);
  
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(y(i)) || isnan(y2(i)))
      {
        cout << "Mlt incorrect" << endl;
        abort();
      }
  
  // checking MltAdd
  SetComplexReal(to_num<Treal>("2.2"), coef_mlt);
  T coef_mlt_beta;
  SetComplexReal(to_num<Treal>("1.2"), coef_mlt_beta);
  y.FillRand(); y *= 1e-9; y2 = y;
  MltAdd(coef_mlt, A, x, coef_mlt_beta, y);
  MltAdd(coef_mlt, B, x, coef_mlt_beta, y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(y(i)) || isnan(y2(i)))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }
  
  // checking GetLU/ SolveLU
  Mlt(A, x, y);
  GetLU(A, Cband, true);
  SolveLU(Cband, y);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - x(i)) > threshold) || isnan(x(i)) || isnan(y(i)))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }
  
  // checking GetLU/SolveLU of Lapack interface
  //n = 26;
  IVect pivot(n);
  Mlt(A, x, y);
  GetLU(A, pivot);
  SolveLU(A, pivot, y);

  for (int i = 0; i < n; i++)
    if ((abs(y(i) - x(i)) > threshold) || isnan(x(i)) || isnan(y(i)))
      {
	DISP(x(i)); DISP(y(i)); DISP(abs(y(i)-x(i)));
        cout << "Lapack's GetLU/SolveLU incorrect" << endl;
        abort();
      }
}

template<class T>
void CheckArrowMatrix(Matrix<T, General, ArrowCol>& A, int n, int kl, int ku, int nrow, int ncol)
{
  A.Reallocate(n, n, kl, ku, nrow, ncol);
  
  // checking GetM/GetN
  if ( (A.GetM() != n) || (A.GetN() != n) )
    {
      cout << "GetM/GetN incorrect" << endl;
      abort();
    }
  
  // checking GetKL, GetKU
  if ( (A.GetKL() != kl) || (A.GetKU() != ku) )
    {
      cout << "GetKL/GetKU incorrect" << endl;
      abort();
    }

  // checking GetNbLastRow, GetNbLastCol
  if ( (A.GetNbLastRow() != nrow) || (A.GetNbLastCol() != ncol) )
    {
      cout << "GetNbLastRow/GetNbLastCol incorrect" << endl;
      abort();
    }
  
  // checking Clear
  A.Clear();
  
  if ( (A.GetM() != 0) || (A.GetN() != 0) )
    {
      cout << "Clear incorrect" << endl;
      abort();
    }
  
  // checking Zero
  A.Reallocate(n, n, kl, ku, nrow, ncol);

  T zero; SetComplexZero(zero);
  T one; SetComplexOne(one);
  typedef typename ClassComplexType<T>::Treal Treal;
  
  A.Zero();
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j)) > threshold) || isnan(A(i, j)))
        {
          cout << "Zero incorrect" << endl;
          abort();
        }
  
  // checking Get/Val
  Matrix<T, General, RowMajor> B(n, n);
  B.FillRand(); Mlt(Treal(1e-9), B);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	  {
	    A.Get(i, j) = B(i, j);
	  }
	else
	  B(i, j) = zero;
      }
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (A(i, j) != B(i, j))
        {
          cout << "Get incorrect" << endl;
          abort();
        }
  
  A.Fill(zero);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	A.Val(i, j) = B(i, j);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
        {
          cout << "Val incorrect" << endl;
          abort();
        }
  
  // checking AddInteraction
  Matrix<T, General, RowMajor> C(n, n);
  C.FillRand(); Mlt(Treal(1e-9), C);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	  A.AddInteraction(i, j, C(i, j));
	else
	  C(i, j) = zero;
      }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j) - C(i, j)) > threshold) || isnan(A(i, j)) || isnan(B(i, j)) || isnan(C(i, j)))
        {
          cout << "AddInteraction incorrect" << endl;
          abort();
        }

  // checking AddInteractionRow
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
	if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	  A.Val(i, j) = B(i, j);
  
  IVect num(n); Vector<T> val(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	int nb_val = 0;
	if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	  {
	    num(nb_val) = j;
	    val(nb_val) = C(i, j);
	    nb_val++;
	  }
	
	A.AddInteractionRow(i, nb_val, num, val);
      }

  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j) - C(i, j)) > threshold) || isnan(A(i, j)) || isnan(B(i, j)) || isnan(C(i, j)))
        {
          cout << "AddInteractionRow incorrect" << endl;
          abort();
        }
  
  // checking ClearRow and operator *=
  T coef_alpha;
  SetComplexReal(to_num<Treal>("1.3"), coef_alpha);
  A.ClearRow(5);
  A *= coef_alpha;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	T val_ref = coef_alpha*(B(i, j) + C(i, j));
	if (i == 5)
	  val_ref = zero;
	
	if ((abs(A(i, j) - val_ref) > threshold) || isnan(A(i, j)) || isnan(val_ref))
	  {
	    cout << "ClearRow incorrect" << endl;
	    abort();
	  }
      }
  
  // checking SetIdentity
  A.SetIdentity();
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        T val_ex = zero;
        if (i == j)
          val_ex = one;
        
        if ((abs(A(i, j)-val_ex) > threshold) || isnan(A(i, j)) || isnan(val_ex))
          {
	    DISP(i); DISP(j); DISP(A(i, j));
	    A.WriteText("test.dat");
            cout << "SetIdentity incorrect" << endl;
            abort();
          }
      }
  
  // checking Fill
  SetComplexReal(to_num<Treal>("2.5"), coef_alpha);
  A.Fill(coef_alpha);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	if ((abs(A(i, j) - coef_alpha) > threshold) || isnan(A(i, j)) || isnan(coef_alpha))
	  {
	    cout << "Fill incorrect" << endl;
	    abort();
	  }
  
  // checking FillRand, Write, WriteText
  A.FillRand(); A *= Treal(1e-9);
  
  A.Write("band_bin.dat");
  
  A.WriteText("band.dat");
  
  // checking conversion from sparse matrices to band matrices
  /*Matrix<T, General, RowSparse> Acsr(n, n);
  Matrix<T, General, ArrayRowSparse> Asp;
  Acsr.FillRand(80); Mlt(1e-9, Acsr);
  Acsr.WriteText("mat_csr.dat");  
  Copy(Acsr, Asp);
  Copy(Asp, A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (Acsr(i, j) != A(i, j))
	{
	  cout << "Copy incorrect" << endl;
	  abort();
	}
  */
  
  // checking Add
  Matrix<T, General, ArrowCol> Cband;
  Cband.Reallocate(n, n, kl, ku, nrow, ncol);
  A.Reallocate(n, n, kl, ku, nrow, ncol);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	{
	  A.Get(i, j) = B(i, j);
	  Cband.Get(i, j) = C(i, j);
	}

  SetComplexReal(to_num<Treal>("0.9"), coef_alpha);
  Add(coef_alpha, Cband, A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j) - coef_alpha*C(i, j)) > threshold) || isnan(A(i, j)) || isnan(C(i, j)) || isnan(B(i, j)))
	{
	  cout << "Add incorrect" << endl;
	  abort();
	}
  
  // checking Mlt
  A.Reallocate(n, n, kl, ku, nrow, ncol);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ( ((i <= j+kl) && (i >=j-ku)) || (i >= n-nrow) || (j >= n-ncol))
	A.Get(i, j) = B(i, j);
  
  SetComplexReal(to_num<Treal>("3.2"), coef_alpha);
  Mlt(coef_alpha, A);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - coef_alpha*B(i, j)) > threshold) || isnan(A(i, j)) || isnan(B(i, j)))
	{
	  cout << "Mlt incorrect" << endl;
	  abort();
	}
  
  Mlt(one/coef_alpha, A);
  Vector<T> x(n), y(n), y2(n);
  x.FillRand(); x *= 1e-9;
  Mlt(A, x, y);
  Mlt(B, x, y2);
  
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(y(i)) || isnan(y2(i)))
      {
        cout << "Mlt incorrect" << endl;
        abort();
      }
  
  // checking MltAdd
  SetComplexReal(to_num<Treal>("2.2"), coef_alpha);
  T coef_beta; SetComplexReal(to_num<Treal>("1.2"), coef_beta);
  y.FillRand(); y *= 1e-9; y2 = y;
  MltAdd(coef_alpha, A, x, coef_beta, y);
  MltAdd(coef_alpha, B, x, coef_beta, y2);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - y2(i)) > threshold) || isnan(y(i)) || isnan(y2(i)))
      {
        cout << "MltAdd incorrect" << endl;
        abort();
      }
  
  // checking GetLU/ SolveLU
  Mlt(A, x, y);
  GetLU(A, Cband, true);
  SolveLU(Cband, y);
  for (int i = 0; i < n; i++)
    if ((abs(y(i) - x(i)) > threshold) || isnan(y(i)) || isnan(x(i)))
      {
        cout << "GetLU/SolveLU incorrect" << endl;
        abort();
      }
}


int main(int argc, char** argv)
{
  bool overall_success = true;
  
  cout.precision(16);
  
  threshold = 1e-11;
  
  {
    // testing band-matrices
    Matrix<Real_wp, General, BandedCol> A;
    int n = 20, kl = 2, ku = 3;
    
    CheckBandMatrix(A, n, kl, ku);
  }
  
  {
    // testing band-matrices
    Matrix<Real_wp, General, BandedCol> A;
    int n = 25, kl = 4, ku = 2;
    
    CheckBandMatrix(A, n, kl, ku);
  }

  {
    // testing arrow matrices
    Matrix<Real_wp, General, ArrowCol> A;
    int n = 22, kl = 3, ku = 4, nb_last_row = 5, nb_last_col = 3;
    
    CheckArrowMatrix(A, n, kl, ku, nb_last_row, nb_last_col);
  }

  {
    // testing arrow matrices
    Matrix<Real_wp, General, ArrowCol> A;
    int n = 27, kl = 5, ku = 3, nb_last_row = 2, nb_last_col = 5;
    
    CheckArrowMatrix(A, n, kl, ku, nb_last_row, nb_last_col);
  }


  {
    // testing band-matrices
    Matrix<Complex_wp, General, BandedCol> A;
    int n = 20, kl = 2, ku = 3;
    
    CheckBandMatrix(A, n, kl, ku);
  }

  {
    // testing band-matrices
    Matrix<Complex_wp, General, BandedCol> A;
    int n = 25, kl = 4, ku = 2;
    
    CheckBandMatrix(A, n, kl, ku);
  }

  {
    // testing arrow matrices
    Matrix<Complex_wp, General, ArrowCol> A;
    int n = 22, kl = 3, ku = 4, nb_last_row = 5, nb_last_col = 3;
    
    CheckArrowMatrix(A, n, kl, ku, nb_last_row, nb_last_col);
  }

  {
    // testing arrow matrices
    Matrix<Complex_wp, General, ArrowCol> A;
    int n = 27, kl = 5, ku = 3, nb_last_row = 2, nb_last_col = 5;
    
    CheckArrowMatrix(A, n, kl, ku, nb_last_row, nb_last_col);
  }
  
  if (overall_success)
    cout << "All tests passed successfully" << endl;
  
  return 0;
}
