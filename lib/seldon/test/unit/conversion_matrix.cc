#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<double> VectReal_wp;
typedef Vector<complex<double> > VectComplex_wp;

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
int largeur = 50;

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
void GenerateRandomVector(Vector<T>& x, int m)
{
  x.Reallocate(m);
  for (int i = 0; i < m; i++)
    GetRand(x(i));  
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRand(x);
      A.Set(i, j, x);
      if (IsSymmetricMatrix(A))
        A.Set(j, i, x);
    }
}

template<class T, class Prop1, class Storage1,
         class T2, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop1, Storage1>& A,
		 const Matrix<T2, Prop2, Storage2>& B)
{
  if ( (A.GetM() == 0) || (A.GetN() == 0)  || (A.GetDataSize() == 0) )
    return false;
  
  if ( (A.GetM() != B.GetM()) || (A.GetN() != B.GetN()) )
    return false;
  
  typename Matrix<T, Prop1, Storage1>::entry_type val, valB;
  for (int i = 0; i < A.GetM(); i++)
    for (int j = max(0, i-largeur); j < min(A.GetN(), i+largeur); j++)
      {
        val = A(i, j);
        valB = B(i, j);
	Real_wp normA = 1;
	if (abs(val) > 1)
	  normA = abs(val);
	
        if ( (abs(val - valB) > threshold*normA ) || isnan(val) || isnan(valB) || isnan(normA))
          {
            DISP(A(i, j)); DISP(B(i, j)); DISP(i); DISP(j);
            return false;
          }
      }
  
  return true;
}


template<class T, class Prop, class Storage, class Allocator>
void CheckRealUnsymmetric(Matrix<T, Prop, Storage, Allocator>& A)
{  
  Matrix<T, Prop, Storage, Allocator> B;
  
  // testing if matrices are correctly read
  ReadHarwellBoeing("test/matrix/market/fidap014.rua", A);
    
  B.ReadText("test/matrix/market/fidap014.dat");
  if ((B.GetM() == 0) || (A.GetM() == 0) || (A.GetDataSize() == 0) || (B.GetDataSize() == 0))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, B))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }        
  
  B.Clear();
  ReadMatrixMarket("test/matrix/market/fidap014.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "ReadMatrixMarket incorrect" << endl;
      abort();
    }
  
  // then testing if matrices are correctly written
  B.WriteText("toto.dat");
  WriteMatrixMarket(B, "toto.mtx");
  WriteHarwellBoeing(B, "toto.rua");
  
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "WriteText incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("toto.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteMatrixMarket incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadHarwellBoeing("toto.rua", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteHarwellBoeing incorrect" << endl;
      abort();
    }

  // testing rectangular matrix
  A.Clear(); B.Clear();
  ReadHarwellBoeing("test/matrix/market/well1033.rra", A);
  
  B.ReadText("test/matrix/market/well1033.dat");
  if ((B.GetM() == 0) || (A.GetM() == 0) || (A.GetDataSize() == 0) || (B.GetDataSize() == 0))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, B))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }        
  
  B.Clear();
  ReadMatrixMarket("test/matrix/market/well1033.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "ReadMatrixMarket incorrect" << endl;
      abort();
    }
  
  // then testing if matrices are correctly written
  B.WriteText("toto.dat");
  WriteMatrixMarket(B, "toto.mtx");
  WriteHarwellBoeing(B, "toto.rua");
  
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "WriteText incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("toto.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteMatrixMarket incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadHarwellBoeing("toto.rua", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteHarwellBoeing incorrect" << endl;
      abort();
    }  
}

template<class T, class Prop, class Storage, class Allocator>
void CheckComplexUnsymmetric(Matrix<T, Prop, Storage, Allocator>& A)
{  
  Matrix<T, Prop, Storage, Allocator> B;
  ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
  
  B.ReadText("test/matrix/market/mhd1280a.dat", true);
  if ((B.GetM() == 0) || (A.GetM() == 0) || (B.GetDataSize() == 0) || (A.GetDataSize() == 0))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, B))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("test/matrix/market/mhd1280a.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "ReadMatrixMarket incorrect" << endl;
      abort();
    }
  
  // then testing if matrices are correctly written
  B.WriteText("toto.dat");
  WriteMatrixMarket(B, "toto.mtx");
  WriteHarwellBoeing(B, "toto.rua");
  
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "WriteText incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("toto.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteMatrixMarket incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadHarwellBoeing("toto.rua", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteHarwellBoeing incorrect" << endl;
      abort();
    }
}
  

template<class T, class Prop, class Storage, class Allocator>
void CheckRealSymmetric(Matrix<T, Prop, Storage, Allocator>& A)
{
  Matrix<T, Prop, Storage, Allocator> B;
  // testing if matrices are correctly read
  ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
  
  B.ReadText("test/matrix/market/bcsstk14.dat");
  if ((B.GetM() == 0) || (A.GetM() == 0) || (A.GetDataSize() == 0) || (B.GetDataSize() == 0))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, B))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }        
  
  B.Clear();
  ReadMatrixMarket("test/matrix/market/bcsstk14.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "ReadMatrixMarket incorrect" << endl;
      abort();
    }
  
  // then testing if matrices are correctly written
  B.WriteText("toto.dat");
  WriteMatrixMarket(B, "toto.mtx");
  WriteHarwellBoeing(B, "toto.rua");
  
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "WriteText incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("toto.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteMatrixMarket incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadHarwellBoeing("toto.rua", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteHarwellBoeing incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckComplexSymmetric(Matrix<T, Prop, Storage, Allocator>& A)
{
  Matrix<T, Prop, Storage, Allocator> B;
  ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
  
  B.ReadText("test/matrix/market/dwg961a.dat", true);
  if ((B.GetM() == 0) || (A.GetM() == 0) || (B.GetDataSize() == 0) || (A.GetDataSize() == 0))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  if (!EqualMatrix(A, B))
    {
      cout << "ReadHarwellBoeing incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("test/matrix/market/dwg961a.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "ReadMatrixMarket incorrect" << endl;
      abort();
    }
  
  // then testing if matrices are correctly written
  B.WriteText("toto.dat");
  WriteMatrixMarket(B, "toto.mtx");
  WriteHarwellBoeing(B, "toto.rua");
  
  B.Clear();
  B.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "WriteText incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadMatrixMarket("toto.mtx", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteMatrixMarket incorrect" << endl;
      abort();
    }
  
  B.Clear();
  ReadHarwellBoeing("toto.rua", B);
  if (!EqualMatrix(A, B))
    {
      cout << "WriteHarwellBoeing incorrect" << endl;
      abort();
    } 
}

template<class T, class Prop, class Storage, class Allocator>
void CheckSymmetryPattern(const Matrix<T, Prop, Storage, Allocator>& A)
{
  /*T zero, val;
  SetComplexZero(zero);
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if (A(i, j) != zero)
	{
	  val = A.Val(j, i);
          if (val != A.Val(i, j))
            {
              cout << "Pattern not symmetric" << endl;
              abort();
            }
	}
  */
}

int main(int argc, char** argv)
{
  
  cout.precision(15);
  threshold = 1e-11;
  //srand(time(NULL));
  //srand(0);

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, General, ColSparse> B;
    Matrix<Real_wp, General, ColSparse> A;
    ReadHarwellBoeing("test/matrix/market/fidap014.rua", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }
#endif

  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, General, RowSparse> B;
    Matrix<Real_wp, General, RowSparse> A;
    ReadHarwellBoeing("test/matrix/market/fidap014.rua", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, General, ArrayRowSparse> B;
    Matrix<Real_wp, General, ArrayRowSparse> A;
    ReadHarwellBoeing("test/matrix/market/fidap014.rua", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, General, ArrayColSparse> B;
    Matrix<Real_wp, General, ArrayColSparse> A;
    ReadHarwellBoeing("test/matrix/market/fidap014.rua", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }
#endif

  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> B;
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, Symmetric, RowSymSparse> B;
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> B;
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }

  {
    // testing conversion to coordinate matrix
    Matrix<Real_wp, Symmetric, ColSymSparse> B;
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
    
    IVect IndRow, IndCol; VectReal_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    
  
  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, General, ColComplexSparse> B;
    Matrix<Complex_wp, General, ColComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }
#endif

  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, General, RowComplexSparse> B;
    Matrix<Complex_wp, General, RowComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, General, ArrayColComplexSparse> B;
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    
#endif
  
  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, General, ArrayRowComplexSparse> B;
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> B;
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
    
    B.Reallocate(A.GetM(), A.GetN());
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    
#endif
  
  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> B;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
    
    B.Reallocate(A.GetM(), A.GetN());
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> B;
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
    
    B.Reallocate(A.GetM(), A.GetN());
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }
#endif

  {
    // testing conversion to coordinate matrix
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> B;
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
    
    B.Reallocate(A.GetM(), A.GetN());
    IVect IndRow, IndCol; VectComplex_wp Val;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    ConvertMatrix_from_Coordinates(IndRow, IndCol, Val, B);
    if (!EqualMatrix(A, B))
      {
	cout << "ConvertMatrix_to_Coordinates incorrect" << endl;
	abort();
      }
  }    
  
  {
    // testing ConvertToCSC for unsymmetric matrices
    int m = 52, n = 52, nnz = 400;
    Matrix<Real_wp, General, RowSparse> A;
    Matrix<Real_wp, General, ColSparse> B, C;
    GenerateRandomMatrix(A, m, n, nnz);
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectReal_wp Val;
    General sym;
    ConvertToCSC(A, sym, Ptr, IndRow, Val, true);
    
    B.SetData(m, n, Val, Ptr, IndRow);
    CheckSymmetryPattern(B);
    
    Copy(A, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
    
    Matrix<Real_wp, General, ArrayRowSparse> A2;
    Copy(A, A2);
    
    B.Clear(); C.Clear();
    ConvertToCSC(A2, sym, Ptr, IndRow, Val, true);
    
    B.SetData(m, n, Val, Ptr, IndRow);
    CheckSymmetryPattern(B);
    
    Copy(A2, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

#ifndef SELDON_WITH_COMPILED_LIBRARY
    Matrix<Real_wp, General, ColSparse> A3;
    Copy(A, A3);
    
    B.Clear(); C.Clear();
    ConvertToCSC(A3, sym, Ptr, IndRow, Val, true);

    B.SetData(m, n, Val, Ptr, IndRow);
    CheckSymmetryPattern(B);
    
    Copy(A3, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, General, ArrayColSparse> A4;
    Copy(A, A4);
    
    B.Clear(); C.Clear();
    ConvertToCSC(A4, sym, Ptr, IndRow, Val, true);

    B.SetData(m, n, Val, Ptr, IndRow);
    CheckSymmetryPattern(B);
    
    Copy(A4, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
#endif
  } 

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing ConvertToCSC for symmetric matrices
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    Matrix<Real_wp, General, ColSparse> B, B2;
    Matrix<Real_wp, Symmetric, ColSymSparse> C, C2;
    
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
    int m = A.GetM(), n = m;
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectReal_wp Val;
    General unsym; Symmetric sym;
    ConvertToCSC(A, unsym, Ptr, IndRow, Val);
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A, sym, Ptr, IndRow, Val);
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A, C2);     Copy(A, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
    
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A2;
    Copy(A, A2);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSC(A2, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A2, sym, Ptr, IndRow, Val);
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A2, C2);     Copy(A2, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, Symmetric, ColSymSparse> A3;
    Copy(A, A3);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSC(A3, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A3, sym, Ptr, IndRow, Val);
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A3, C2);     Copy(A3, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A4;
    Copy(A, A4);

    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSC(A4, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A4, sym, Ptr, IndRow, Val);
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A4, C2);     Copy(A4, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
  } 
#endif
  
  {
    // testing ConvertToCSR for unsymmetric matrices
    int m = 52, n = 55, nnz = 400;
    Matrix<Real_wp, General, RowSparse> A;
    Matrix<Real_wp, General, RowSparse> B, C;
    GenerateRandomMatrix(A, m, n, nnz);
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectReal_wp Val;
    General sym;
    ConvertToCSR(A, sym, Ptr, IndRow, Val);
    
    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
    
    Matrix<Real_wp, General, ArrayRowSparse> A2;
    Copy(A, A2);
    
    B.Clear(); C.Clear();
    ConvertToCSR(A2, sym, Ptr, IndRow, Val);
    
    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A2, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

#ifndef SELDON_WITH_COMPILED_LIBRARY
    Matrix<Real_wp, General, ColSparse> A3;
    Copy(A, A3);
    
    B.Clear(); C.Clear();
    ConvertToCSR(A3, sym, Ptr, IndRow, Val);

    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A3, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, General, ArrayColSparse> A4;
    Copy(A, A4);
    
    B.Clear(); C.Clear();
    ConvertToCSR(A4, sym, Ptr, IndRow, Val);

    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A4, C);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
#endif
  } 

  {
    // testing ConvertToCSR for symmetric matrices
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    Matrix<Real_wp, General, RowSparse> B, B2;
    Matrix<Real_wp, Symmetric, RowSymSparse> C, C2;
    
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", A);
    int m = A.GetM(), n = m;
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectReal_wp Val;
    General unsym; Symmetric sym;
    ConvertToCSR(A, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A, sym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A, B2); Copy(A, C2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
    
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A2;
    Copy(A, A2);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSR(A2, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A2, sym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A2, B2); Copy(A2, C2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

#ifndef SELDON_WITH_COMPILED_LIBRARY
    Matrix<Real_wp, Symmetric, ColSymSparse> A3;
    Copy(A, A3);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSR(A3, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A3, sym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A3, B2); Copy(A3, C2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A4;
    Copy(A, A4);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSR(A4, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A4, sym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A4, B2); Copy(A4, C2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, C) ||
        !EqualMatrix(A, B2) || !EqualMatrix(A, C2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
#endif
  } 

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing other unsymmetric conversions (function Copy)
    int m = 52, n = 48, nnz = 400;
    Matrix<Real_wp, General, RowSparse> A;
    Matrix<Real_wp, General, ColSparse> B;
    GenerateRandomMatrix(A, m, n, nnz);
    
    
    Matrix<Real_wp, General, ArrayRowSparse> A2;
    Copy(A, A2); Copy(A, B);
    if (!EqualMatrix(A, A2) || !EqualMatrix(A, B))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, General, ArrayColSparse> B2;
    Copy(B, B2);
    if (!EqualMatrix(B, B2))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    A2.Clear();
    Copy(B, A2);
    if (!EqualMatrix(B, A2))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    A2.Clear();
    Copy(B2, A2);
    if (!EqualMatrix(B2, A2))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }    

    B2.Clear();
    Copy(A2, B2);
    if (!EqualMatrix(B2, A2))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }    
  }
#endif
  
  {
    // testing other symmetric conversions (function Copy)
    int m = 51, n = m, nnz = 400;
    Matrix<Real_wp, General, RowSparse> A, B;
    Matrix<Real_wp, General, ArrayRowSparse> Ar, Br;
    Matrix<Real_wp, Symmetric, RowSymSparse> As;
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> C, Bs;
    GenerateRandomMatrix(A, m, n, nnz);
    
    Real_wp one(1);
    B = A; 
    Transpose(B);
    Copy(A, Ar); Copy(B, Br);
    Add(one, Ar, Br);
    
    Copy(Br, C);
    Copy(C, As); Copy(As, Bs);
    if (!EqualMatrix(Br, C) || !EqualMatrix(As, Bs) || !EqualMatrix(As, C))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Bs.Clear();
#ifndef SELDON_WITH_COMPILED_LIBRARY    
    Matrix<Real_wp, Symmetric, ColSymSparse> As2;
    Copy(As, As2); Copy(As2, Bs);
    if (!EqualMatrix(As, As2) || !EqualMatrix(As2, Bs))
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> Bs2;
    Copy(As, Bs2);
    if (!EqualMatrix(As, Bs2) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Bs.Clear();
    Copy(Bs2, Bs);
    if (!EqualMatrix(Bs, Bs2) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
        
    Ar.Clear();    
    Copy(Bs, Ar);
    if (!EqualMatrix(Bs, Ar) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Real_wp, General, ArrayColSparse> Br2;
    Copy(Bs, Br2);
    if (!EqualMatrix(Bs, Br2) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
#endif
  }

  {
    // testing GetSymmetricPattern
    int m = 49, n = m, nnz = 400;
    Matrix<int, General, RowSparse> A;
    A.Reallocate(m, n);
    A.FillRand(nnz);
    
    Matrix<int, Symmetric, RowSymSparse> B;
    GetSymmetricPattern(A, B);
    
    A.Fill(1);
    B.Fill(1);
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
	{
	  if (A(i, j) != 0)
	    {
	      if ( (B(i, j) != 1) || (B(j, i) != 1) )
		{
		  cout << "GetSymmetricPattern incorrect" << endl;
		  abort();
		}
	    }
	  
	  if (B(i, j) == 1)
	    {
	      if ( (A(i, j) == 0) && (A(j, i) == 0) )
		{
		  cout << "GetSymmetricPattern incorrect" << endl;
		  abort();
		}
	    }
	}
  }
  
  {
    // testing conversion Dense-Sparse matrices
    Matrix<Real_wp, General, RowSparse> A;
    Matrix<Real_wp, General, ArrayRowSparse> A2;
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> As;
    Matrix<Real_wp, General, RowMajor> B;
    Matrix<Real_wp, Symmetric, RowSymPacked> Bs;
    
    ReadHarwellBoeing("test/matrix/market/fidap014.rua", A);
    Copy(A, A2);
    Copy(A, B);
    if (!EqualMatrix(A, B) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    B.Clear(); Copy(A2, B);
    if (!EqualMatrix(A2, B) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
        
    ReadHarwellBoeing("test/matrix/market/bcsstk14.rsa", As);
    Copy(As, Bs);
    if (!EqualMatrix(As, Bs) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Matrix<Real_wp, Symmetric, RowSymSparse> As2;
    ConvertToSparse(Bs, As2, threshold);
    if (!EqualMatrix(As2, Bs) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    A.Clear();
    ConvertToSparse(B, A, threshold);
    ConvertToSparse(B, A2, threshold);
    if (!EqualMatrix(A2, B) || !EqualMatrix(A, B) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }    
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    // testing ConvertToCSC for complex matrices
    Matrix<Complex_wp, General, RowComplexSparse> A;
    Matrix<Complex_wp, General, ColSparse> B, B2;
    
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    int m = A.GetM(), n = A.GetN();

    Vector<long> Ptr; Vector<int> IndRow;
    VectComplex_wp Val;
    General unsym;
    ConvertToCSC(A, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
    
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A2;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A2);
    
    B.Clear(); B2.Clear();
    ConvertToCSC(A2, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    Copy(A2, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, General, ColComplexSparse> A3;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A3);
    
    B.Clear(); B2.Clear();
    ConvertToCSC(A3, unsym, Ptr, IndRow, Val);
    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A3, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, General, ArrayColComplexSparse> A4;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A4);
    
    B.Clear(); B2.Clear();
    ConvertToCSC(A4, unsym, Ptr, IndRow, Val);
    B.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A4, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
  } 

  {
    // testing ConvertToCSC for symmetric complex matrices
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    Matrix<Complex_wp, Symmetric, ColSymSparse> B, B2;
    Matrix<Complex_wp, General, ColSparse> C, C2;
    
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
    int m = A.GetM(), n = m;
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectComplex_wp Val;
    Symmetric sym; General unsym;
    ConvertToCSC(A, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A, B2); Copy(A, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
    
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A2;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A2);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSC(A2, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A2, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A2, B2); Copy(A2, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A3;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A3);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSC(A3, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A3, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A3, B2); Copy(A3, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A4;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A4);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSC(A4, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSC(A4, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A4, B2); Copy(A4, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSC incorrect" << endl;
	abort();
      }
  } 
#endif
  
  {
    // testing ConvertToCSR for unsymmetric matrices
    Matrix<Complex_wp, General, RowComplexSparse> A;
    Matrix<Complex_wp, General, RowSparse> B, B2;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    int m = A.GetM(), n = A.GetN();
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectComplex_wp Val;
    General unsym;
    ConvertToCSR(A, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);    
    
    Copy(A, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
    
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A2;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A2);
    
    B.Clear(); B2.Clear();
    ConvertToCSR(A2, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);    
    
    Copy(A2, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

#ifndef SELDON_WITH_COMPILED_LIBRARY
    Matrix<Complex_wp, General, ColComplexSparse> A3;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A3);
    
    B.Clear(); B2.Clear();
    ConvertToCSR(A3, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);    
    
    Copy(A3, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, General, ArrayColComplexSparse> A4;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A4);
    
    B.Clear(); B2.Clear();
    ConvertToCSR(A4, unsym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);    
    
    Copy(A4, B2);
    if (!EqualMatrix(A, B) || !EqualMatrix(A, B2) )
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
#endif
  } 

  {
    // testing ConvertToCSR for symmetric matrices
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    Matrix<Complex_wp, Symmetric, RowSymSparse> B, B2;
    Matrix<Complex_wp, General, RowSparse> C, C2;
        
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A);
    int m = A.GetM(), n = m;
    
    Vector<long> Ptr; Vector<int> IndRow;
    VectComplex_wp Val;
    Symmetric sym; General unsym;
    ConvertToCSR(A, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A, B2); Copy(A, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
        
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A2;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A2);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSR(A2, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A2, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A2, B2); Copy(A2, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

#ifndef SELDON_WITH_COMPILED_LIBRARY
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A3;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A3);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSR(A3, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A3, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A3, B2); Copy(A3, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A4;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", A4);
    
    B.Clear(); C.Clear(); B2.Clear(); C2.Clear();
    ConvertToCSR(A4, sym, Ptr, IndRow, Val);    
    B.SetData(m, n, Val, Ptr, IndRow);

    ConvertToCSR(A4, unsym, Ptr, IndRow, Val);    
    C.SetData(m, n, Val, Ptr, IndRow);
    
    Copy(A4, B2); Copy(A4, C2);
    if (!EqualMatrix(A, C) || !EqualMatrix(A, B) ||
        !EqualMatrix(A, C2) || !EqualMatrix(A, B2))
      {
	cout << "ConvertToCSR incorrect" << endl;
	abort();
      }
    #endif
  } 

  {
    // testing other conversions (function Copy)
    int m = 65, nnz = 600;
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    Matrix<Complex_wp, General, ArrayRowSparse> B;
    GenerateRandomMatrix(A, m, m, nnz);
    
    Copy(A, B);
    if (!EqualMatrix(A, B))
      {
        cout << "Copy incorrect" << endl;
        abort();
      }
  }
  
  {
    // testing other conversions with complex matrices (function Copy)
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    Matrix<Complex_wp, General, RowComplexSparse> B;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", A);
    
    Copy(A, B);
    if (!EqualMatrix(A, B) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, General, ArrayRowSparse> Ba;    
    Copy(A, Ba);
    if (!EqualMatrix(A, Ba) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> As;
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> Bs;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", As);
    
    Copy(As, Bs);
    if (!EqualMatrix(As, Bs) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Ba.Clear();
    Copy(As, Ba);
    if (!EqualMatrix(As, Ba) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> Ba2;
    
    Copy(As, Ba2);
    if (!EqualMatrix(As, Ba2) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> As2;
    ReadHarwellBoeing("test/matrix/market/dwg961a.csa", As2);
    
    Ba.Clear();
    Copy(As2, Ba);
    if (!EqualMatrix(As2, Ba) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }
    
    Ba2.Clear();
    Copy(As2, Ba2);
    if (!EqualMatrix(As2, Ba2) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

    Matrix<Complex_wp, General, RowComplexSparse> Aunsym;
    ReadHarwellBoeing("test/matrix/market/mhd1280a.cua", Aunsym);
    
    Ba.Clear();
    Copy(Aunsym, Ba);
    if (!EqualMatrix(Aunsym, Ba) )
      {
	cout << "Copy incorrect" << endl;
	abort();
      }

  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing real ColSparse..." << endl;
    Matrix<Real_wp, General, ColSparse> A;
    CheckRealUnsymmetric(A);
  }
    
  {
    cout << "Testing complex ColSparse..." << endl;
    Matrix<Complex_wp, General, ColSparse> A;
    CheckComplexUnsymmetric(A);
  }
#endif

  {
    cout << "Testing real RowSparse..." << endl;
    Matrix<Real_wp, General, RowSparse> A;
    CheckRealUnsymmetric(A);
  }
  
  {
    cout << "Testing complex RowSparse..." << endl;
    Matrix<Complex_wp, General, RowSparse> A;
    CheckComplexUnsymmetric(A);
  }

  {
    cout << "Testing real RowSymSparse..." << endl;
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    CheckRealSymmetric(A);
  }

  {
    cout << "Testing complex RowSymSparse..." << endl;
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    CheckComplexSymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing real ColSymSparse..." << endl;
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    CheckRealSymmetric(A);
  }

  {
    cout << "Testing complex ColSymSparse..." << endl;
    Matrix<Complex_wp, Symmetric, ColSymSparse> A;
    CheckComplexSymmetric(A);
  }
#endif

  {
    cout << "Testing real ArrayRowSparse..." << endl;
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckRealUnsymmetric(A);
  }

  {
    cout << "Testing complex ArrayRowSparse..." << endl;
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckComplexUnsymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing real ArrayColSparse..." << endl;
    Matrix<Real_wp, General, ArrayColSparse> A;
    CheckRealUnsymmetric(A);
  }

  {
    cout << "Testing complex ArrayColSparse..." << endl;
    Matrix<Complex_wp, General, ArrayColSparse> A;
    CheckComplexUnsymmetric(A);
  }
#endif
  
  {
    cout << "Testing real ArrayRowSymSparse..." << endl;
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckRealSymmetric(A);
  }

  {
    cout << "Testing complex ArrayRowSymSparse..." << endl;
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckComplexSymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing real ArrayColSymSparse..." << endl;
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    CheckRealSymmetric(A);
  }

  {
    cout << "Testing complex ArrayColSymSparse..." << endl;
    Matrix<Complex_wp, Symmetric, ArrayColSymSparse> A;
    CheckComplexSymmetric(A);
  }
#endif

  {
    cout << "Testing RowComplexSparse..." << endl;
    Matrix<Complex_wp, General, RowComplexSparse> A;
    CheckComplexUnsymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing ColComplexSparse..." << endl;
    Matrix<Complex_wp, General, ColComplexSparse> A;
    CheckComplexUnsymmetric(A);
  }
#endif

  {
    cout << "Testing RowSymComplexSparse..." << endl;
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    CheckComplexSymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing ColSymComplexSparse..." << endl;
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
    CheckComplexSymmetric(A);
  }
#endif

  {
    cout << "Testing ArrayRowComplexSparse..." << endl;
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    CheckComplexUnsymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing ArrayColComplexSparse..." << endl;
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
    CheckComplexUnsymmetric(A);
  }
#endif

  {
    cout << "Testing ArrayRowSymComplexSparse..." << endl;
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    CheckComplexSymmetric(A);
  }

#ifndef SELDON_WITH_COMPILED_LIBRARY
  {
    cout << "Testing ArrayColSymComplexSparse..." << endl;
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
    CheckComplexSymmetric(A);
  }
#endif
  
  std::remove("toto.dat");
  std::remove("toto.mtx");
  std::remove("toto.rua");
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
