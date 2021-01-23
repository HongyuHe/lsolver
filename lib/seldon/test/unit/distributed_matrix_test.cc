#include "SeldonLib.hxx"

using namespace Seldon;

typedef Vector<Real_wp> VectReal_wp;

Real_wp threshold;
int rank_processor, nb_processors;

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

template<class T, class Allocator>
void GenerateRandomMatrix(Matrix<T, General, ColMajor, Allocator>& A,
                          int m, int n)
{
  T x;
  A.Reallocate(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        GetRand(x);
        A.Set(i, j, x);
      }
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


template<class T1, class Prop1, class Storage1, class Allocator1,
         class T2, class Prop2, class Storage2, class Allocator2>
bool EqualMatrix(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
                 const Matrix<T2, Prop2, Storage2, Allocator2>& B)
{
  if ( (A.GetM() != B.GetM())  || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j)) > 100.0*threshold) || isnan(abs(A(i, j) - B(i, j))))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          return false;
        }
  
  return true;
}


void GenerateRandomPermutation(int n, IVect& permut)
{
  Vector<bool> NumUsed(n);
  NumUsed.Fill(false);
  permut.Reallocate(n);
  permut.Fill(-1);
  int nb = 0;
  // premiere iteration
  for (int i = 0; i < n; i++)
    {
      int i2 = rand()%n;
      if (!NumUsed(i2))
        {
          NumUsed(i2) = true;
          permut(i) = i2;
          nb++;
        }
    }
  
  while (nb < n)
    {
      // on recupere les numeros non-selectionnes
      IVect non_selec(n-nb);
      int k = 0;
      for (int i = 0; i < n; i++)
        if (!NumUsed(i))
          non_selec(k++) = i;
      
      // iteration suivante
      for (int i = 0; i < n; i++)
        if (permut(i) == -1)
          {
            int i2 = rand()%(n-nb);
            int j = non_selec(i2);
            if (!NumUsed(j))
              {
                NumUsed(j) = true;
                permut(i) = j;
                nb++;
              }
          }
    }
}

void PickRandomInteger(int nb, IVect& num)
{
  if (nb == num.GetM())
    return;
  
  IVect permut;
  GenerateRandomPermutation(num.GetM(), permut);
  
  IVect old_num(num);
  num.Reallocate(nb);
  for (int i = 0; i < nb; i++)
    num(i) = old_num(permut(i));
}

template<class T>
void PickRandomPartition(const T& x, int n, Vector<T>& eval)
{
  eval.Reallocate(n);
  Real_wp zero, sum; SetComplexZero(zero);
  sum = zero; VectReal_wp coef(n);
  for (int i = 0; i < n; i++)
    {
      GetRand(coef(i));
      sum += coef(i);
    }
  
  Mlt(Real_wp(1)/sum, coef);
  for (int i = 0; i < n; i++)
    eval(i) = coef(i)*x;
}

template<class MatrixSparse, class T>
void AddInteraction(MatrixSparse& A, int i, int j, const T& x,
                    const IVect& Glob_to_local, int proc_row, int proc_col)
{
  int iloc = Glob_to_local(i);
  int jloc = Glob_to_local(j);
  bool sym = IsSymmetricMatrix(A);
  if (i == j)
    sym = false;
  
  if (iloc >= 0)
    {
      if (jloc >= 0)
        {
          A.AddInteraction(iloc, jloc, x);
          if (sym)
            A.AddInteraction(jloc, iloc, x);
        }
      else
        {
          A.AddDistantInteraction(iloc, j, proc_col, x); 
          if (sym)
            A.AddRowDistantInteraction(j, iloc, proc_col, x); 
        }
    }
  else
    {
      A.AddRowDistantInteraction(i, jloc, proc_row, x);
      if (sym)
        A.AddDistantInteraction(jloc, i, proc_row, x);
    }
}

template<class MatrixSparse, class MatrixSparse2>
void CheckMatrixMlt(const MatrixSparse& A,
                    const MatrixSparse2& Aref,
                    const string& fct_name, bool all_col = false)
{
  typedef typename MatrixSparse::entry_type T;
  int m = A.GetGlobalM();
  Vector<T> Xref, Yref;
  GenerateRandomVector(Xref, m);
  GenerateRandomVector(Yref, m);

  T alpha, beta;
  GetRand(alpha);   GetRand(beta);
  
  Vector<T> X(A.GetM()), Y(A.GetM());
  const IVect& global = A.GetGlobalRowNumber();
  for (int i = 0; i < global.GetM(); i++)
    {
      X(i) = Xref(global(i));
      Y(i) = Yref(global(i));
    }
  
  if (all_col)
    {
      Vector<T> Ones(A.GetM()), Ah_Ones(A.GetM());
      for (int j = 0; j < m; j++)
        {
          Ones.Fill(0);
          for (int k = 0; k < global.GetM(); k++)
            if (global(k) == j)
              SetComplexOne(Ones(k));
          
          Mlt(A, Ones, Ah_Ones);
          
          for (int i = 0; i < global.GetM(); i++)
            if ((abs(Ah_Ones(i) - Aref(global(i), j)) > threshold) || isnan(abs(Ah_Ones(i) - Aref(global(i), j))))
              {
                cout << fct_name << " incorrect" << endl;
                DISP(global(i)); DISP(j); DISP(Ah_Ones(i));
                DISP(Aref(global(i), j));
                abort();
              }
        }
    }
  MltAdd(alpha, A, X, beta, Y);
  MltAdd(alpha, Aref, Xref, beta, Yref);
  
  for (int i = 0; i < global.GetM(); i++)
    if ((abs(Yref(global(i)) - Y(i)) > threshold) || isnan(abs(Yref(global(i)) - Y(i))))
      {
        cout << fct_name << " incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(global(i)));
        abort();
      }

  MltAdd(alpha, SeldonTrans, A, X, beta, Y);
  MltAdd(alpha, SeldonTrans, Aref, Xref, beta, Yref);
  
  for (int i = 0; i < global.GetM(); i++)
    if ((abs(Yref(global(i)) - Y(i)) > threshold) || isnan(abs(Yref(global(i)) - Y(i))))
      {
        cout << fct_name << " incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(global(i)));
        abort();
      }


  MltAdd(alpha, SeldonConjTrans, A, X, beta, Y);
  MltAdd(alpha, SeldonConjTrans, Aref, Xref, beta, Yref);
  
  for (int i = 0; i < global.GetM(); i++)
    if ((abs(Yref(global(i)) - Y(i)) > threshold) || isnan(abs(Yref(global(i)) - Y(i))))
      {
        cout << fct_name << " incorrect" << endl;
        DISP(i); DISP(Y(i)); DISP(Yref(global(i)));
        abort();
      }

}

template<class MatrixSeq, class MatrixPar>
void DistributeMatrixProcessor(const MatrixSeq& Aref, MatrixPar& A,
                               const Vector<IVect>& list_proc, const IVect& Glob_to_local)
{
  int n = A.GetM();
  A.Clear(); A.Reallocate(n, n);
  typedef typename MatrixSeq::entry_type T;
  T zero; SetComplexZero(zero);
  for (int i = 0; i < Aref.GetM(); i++)
    {
      int jlow = 0;
      if (IsSymmetricMatrix(Aref))
        jlow = i;
      
      for (int j = jlow; j < Aref.GetM(); j++)
        if (Aref(i, j) != zero)
          {
            IVect num(list_proc(i).GetM() + list_proc(j).GetM());
            for (int k = 0; k < list_proc(i).GetM(); k++)
              num(k) = list_proc(i)(k);
            
            for (int k = 0; k < list_proc(j).GetM(); k++)
              num(list_proc(i).GetM() + k) = list_proc(j)(k);          
            
            RemoveDuplicate(num);
            
            int nb = rand()%num.GetM() + 1;
            PickRandomInteger(nb, num);
            
            IVect num_row(nb), num_col(nb);
            for (int k = 0; k < nb; k++)
              {
                num_row(k) = rand()%list_proc(i).GetM();
                num_col(k) = rand()%list_proc(j).GetM();
              }
            
            Vector<T> eval;
            PickRandomPartition(Aref(i, j), nb, eval);
            
            for (int k = 0; k < nb; k++)
            if (rank_processor == num(k))
              AddInteraction(A, i, j, eval(k), Glob_to_local,
                             list_proc(i)(num_row(k)), list_proc(j)(num_col(k)));
          }
    }
}

template<class T, class Storage, class Allocator, class T1, class Allocator1>
void ScaleLeftMatrix(Matrix<T, Symmetric, Storage, Allocator>&, const Vector<T1, VectFull, Allocator1>&)
{}

template<class T, class Storage, class Allocator, class T1, class Allocator1>
void ScaleRightMatrix(Matrix<T, Symmetric, Storage, Allocator>&, const Vector<T1, VectFull, Allocator1>&)
{}

template<class T, class Prop, class Storage, class Allocator>
void CheckDistributedMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& mat,
                            bool complex_mat = false, bool check_erase = true)
{
  Matrix<T, Prop, Storage, Allocator> Aref;
  typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type T0;
  
  //int m = 4, n = m, nnz = 10;
  int m = 50, n = m, nnz = 300;
  GenerateRandomMatrix(Aref, m, n, nnz);
  T0 x_diag;
  for (int i = 0; i < m; i++)
    {
      GetRand(x_diag);
      Aref.AddInteraction(i, i, x_diag);
    }
  
  // for each dof determining processors that share this dof
  Vector<IVect> list_proc(m);
  for (int i = 0; i < m; i++)
    {
      int nb = rand()%nb_processors+1;
      //int nb = 1;
      
      IVect num(nb);
      for (int j = 0; j < nb; j++)
        num(j) = rand()%nb_processors;
      
      RemoveDuplicate(num);
      nb = num.GetM();
      list_proc(i) = num;
    }
  
  // constructing GlobalRowNumbers
  int nodl = 0, nodl_overlap = 0;
  IVect Glob_to_local(m); Glob_to_local.Fill(-1);
  Vector<int> OverlappedGlobal(m); OverlappedGlobal.Fill(-1);
  IVect NbSharedDofPerProc(nb_processors); NbSharedDofPerProc.Fill(0);
  for (int i = 0; i < list_proc.GetM(); i++)
    for (int j = 0; j < list_proc(i).GetM(); j++)
      if (list_proc(i)(j) == rank_processor)
        {
          Glob_to_local(i) = nodl;
          if (j > 0)
            {
              OverlappedGlobal(i) = list_proc(i)(0);
              nodl_overlap++;
            }
          
          for (int k = 0; k < list_proc(i).GetM(); k++)
            if (j != k)
              NbSharedDofPerProc(list_proc(i)(k))++;
          
          nodl++;
        }
  
  IVect GlobalRowNumbers(nodl);
  for (int i = 0; i < m; i++)
    if (Glob_to_local(i) >= 0)
      GlobalRowNumbers(Glob_to_local(i)) = i;
  
  //DISP(nodl); DISP(GlobalRowNumbers);
  
  // then OverlapRowNumbers/OverlapProcNumbers
  IVect OverlapRowNumbers(nodl_overlap), OverlapProcNumbers(nodl_overlap);
  nodl_overlap = 0;
  for (int i = 0; i < m; i++)
    if (OverlappedGlobal(i) >= 0)
      {
        OverlapRowNumbers(nodl_overlap) = Glob_to_local(i);
        OverlapProcNumbers(nodl_overlap) = OverlappedGlobal(i);
        nodl_overlap++;
      }
  
  //DISP(OverlapRowNumbers); DISP(OverlapProcNumbers);
  // and finally ProcSharingRows, SharingRowNumbers
  int nodl_scalar = nodl; // nb_u = 1
  int nb_proc = 0;
  IVect IndexProc(nb_processors); IndexProc.Fill(-1);
  for (int p = 0; p < nb_processors; p++)
    if (NbSharedDofPerProc(p) > 0)
      {
        IndexProc(p) = nb_proc;
        nb_proc++;
      }
  
  IVect ProcSharingRows(nb_proc);
  Vector<IVect> SharingRowNumbers(nb_proc);
  for (int p = 0; p < nb_processors; p++)
    if (IndexProc(p) >= 0)
      {
        ProcSharingRows(IndexProc(p)) = p;
        SharingRowNumbers(IndexProc(p)).Reallocate(NbSharedDofPerProc(p));
        SharingRowNumbers(IndexProc(p)).Fill(-1);
      }
  
  NbSharedDofPerProc.Fill(0);
  for (int i = 0; i < list_proc.GetM(); i++)
    for (int j = 0; j < list_proc(i).GetM(); j++)
      if (list_proc(i)(j) == rank_processor)
        {
          for (int k = 0; k < list_proc(i).GetM(); k++)
            if (j != k)
              {
                int cpt = NbSharedDofPerProc(list_proc(i)(k));
                int proc = IndexProc(list_proc(i)(k));
                SharingRowNumbers(proc)(cpt) = Glob_to_local(i);
                NbSharedDofPerProc(list_proc(i)(k))++;
              }
        }
  
  for (int i = 0; i < nb_proc; i++)
    {
      //DISP(i); 
      //DISP(ProcSharingRows(i));
      //DISP(SharingRowNumbers(i));
    }

  // distributing Aref to all the processors
  DistributedMatrix<T, Prop, Storage, Allocator> A;
  A.Reallocate(nodl, nodl);
  A.Init(m, &GlobalRowNumbers, &OverlapRowNumbers, &OverlapProcNumbers,
         nodl_scalar, 1, &ProcSharingRows, &SharingRowNumbers, MPI_COMM_WORLD);

  DistributeMatrixProcessor(Aref, A, list_proc, Glob_to_local);
  
  // testing Mlt function
  CheckMatrixMlt(A, Aref, "Mlt", true);

  if (rank_processor == 0)
    Aref.WriteText("mat_ref.dat");
  
  A.WriteText("mat.dat");

  {
    // checking Init when only global row numbers are provided
    DistributedMatrix<T, Prop, Storage, Allocator> Atest;
    IVect GlobalRowNumbersTest, OverlapRowTest, OverlapProcTest;
    IVect ProcSharingTest; Vector<IVect> SharingRowTest;
    GlobalRowNumbersTest = GlobalRowNumbers;
    Atest.Reallocate(nodl, nodl);
    Atest.Init(GlobalRowNumbersTest, OverlapRowTest, OverlapProcTest,
               ProcSharingTest, SharingRowTest, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int i = 0; i < GlobalRowNumbers.GetM(); i++)
      if (GlobalRowNumbers(i) != GlobalRowNumbersTest(i))
        {
          cout << "Incorrect global row number " << endl;
          abort();
        }

    for (int i = 0; i < OverlapRowNumbers.GetM(); i++)
      if (OverlapRowNumbers(i) != OverlapRowTest(i))
        {
          cout << "Incorrect overlap row number " << endl;
          abort();
        }
    
    for (int i = 0; i < OverlapProcTest.GetM(); i++)
      if (OverlapProcNumbers(i) != OverlapProcTest(i))
        {
          cout << "Incorrect overlap proc number " << endl;
          abort();
        }

    for (int i = 0; i < ProcSharingRows.GetM(); i++)
      if (ProcSharingRows(i) != ProcSharingTest(i))
        {
          cout << "Incorrect matching proc number " << endl;
          abort();
        }
    
    for (int i = 0; i < SharingRowNumbers.GetM(); i++)
      for (int j = 0; j < SharingRowNumbers(i).GetM(); j++)
        if (SharingRowNumbers(i)(j) != SharingRowTest(i)(j))
          {
            cout << "Incorrect matching row number " << endl;
            abort();
          }
    
    DistributeMatrixProcessor(Aref, Atest, list_proc, Glob_to_local);
    
    CheckMatrixMlt(Atest, Aref, "Mlt", true);
  }
  
  
  // testing Resize (to be done)

  // testing GetNonZeros/GetDataSize()
  //DISP(Aref.GetNonZeros());
  //DISP(Aref.GetDataSize());
  
  //DISP(A.GetNonZeros());
  //DISP(A.GetDataSize());  
  
  // testing RemoveSmallEntry (to be done)
  
  // testing SetIdentity
  DistributedMatrix<T, Prop, Storage, Allocator> B(A);
  Matrix<T, Prop, Storage, Allocator> Bref(Aref);
  
  CheckMatrixMlt(B, Aref, "Copy constructor");

  B.RemoveSmallEntry(0.1);
  B.SetIdentity();
  Bref.SetIdentity();
  
  CheckMatrixMlt(B, Bref, "SetIdentity");
  
  // testing =
  B = A;
  Bref = Aref;

  if (rank_processor == 0)
    Bref.WriteText("mat_ref_bis.dat");
  
  B.WriteText("mat_bis.dat");
  
  CheckMatrixMlt(B, Bref, "Operator =", true);

  // testing Copy (conversion)
  DistributedMatrix<T0, General, RowSparse> Asparse;
  Copy(A, Asparse);
  
  CheckMatrixMlt(Asparse, Aref, "Copy", true);
  
  // testing Zero
#ifdef MONTJOIE_WITH_MPFR
  T0 zero; SetComplexZero(zero);
  B.Fill(zero);
  Bref.Fill(zero);
#else
  B.Zero();
  Bref.Zero();
#endif
  
  CheckMatrixMlt(B, Bref, "Zero");

  // testing Fill
  B.Fill();
  Bref.Fill();
  
  // CheckMatrixMlt(B, Bref, "Fill");

  T0 alpha, beta; GetRand(alpha); GetRand(beta);
  B.Fill(alpha);
  Bref.Fill(alpha);
  // CheckMatrixMlt(B, Bref, "Fill");
  
  GenerateRandomMatrix(Bref, m, n, nnz+30);
  DistributeMatrixProcessor(Bref, B, list_proc, Glob_to_local);
  
  CheckMatrixMlt(B, Bref, "Test");

  // testing Mlt(scalaire, matrice)
  if (complex_mat)
    {      
      Mlt(realpart(alpha), B);
      Mlt(realpart(alpha), Bref);
    }
  else
    {
      Mlt(alpha, B);
      Mlt(alpha, Bref);
    }

  CheckMatrixMlt(B, Bref, "Mlt", true);
  
  // SOR/GaussSeidel not available
  
  // testing Add
  if (complex_mat)
    {      
      Add(realpart(alpha), Aref, Bref);
      Add(realpart(alpha), A, B);
    }
  else
    {      
      Add(alpha, Aref, Bref);
      Add(alpha, A, B);
    }
    

  if (rank_processor == 0)
    Bref.WriteText("add_ref.dat");
  
  B.WriteText("add.dat");
  
  CheckMatrixMlt(B, Bref, "Add");
  MPI_Barrier(MPI_COMM_WORLD);
  
  // testing MaxAbs
  Real_wp val = MaxAbs(B), val_ref = MaxAbs(Bref);
  
  //DISP(val); DISP(val_ref);
  if (abs(val - val_ref) > threshold)
    {
      //cout << "MaxAbs incorrect" << endl;
      //abort();
    }
  
  // testing GetRowSum
  Vector<Real_wp> sumA, sumAref, sum_colA, sum_colAref;
  
  GetRowSum(sumAref, Aref);
  GetRowSum(sumA, A);
  
  for (int i = 0; i < nodl; i++)
    if ((abs(sumA(i) - sumAref(GlobalRowNumbers(i))) > threshold) || isnan(abs(sumA(i) - sumAref(GlobalRowNumbers(i)))))
      {
        //cout << "GetRowSum incorrect" << endl;
        //DISP(sumA(i)); DISP(sumAref(GlobalRowNumbers(i)));
        //abort();
      }

  GetColSum(sum_colAref, Aref);
  GetColSum(sum_colA, A);
  
  for (int i = 0; i < nodl; i++)
    if ((abs(sum_colA(i) - sum_colAref(GlobalRowNumbers(i))) > threshold) || isnan(abs(sum_colA(i) - sum_colAref(GlobalRowNumbers(i)))))
      {
        //cout << "GetColSum incorrect" << endl;
        //DISP(sum_colA(i)); DISP(sum_colAref(GlobalRowNumbers(i)));
        //abort();
      }

  sumA.Clear(); sum_colA.Clear();
  GetRowColSum(sumA, sum_colA, A);
  
  for (int i = 0; i < nodl; i++)
    if ( (abs(sumA(i) - sumAref(GlobalRowNumbers(i))) > threshold) ||
         (abs(sum_colA(i) - sum_colAref(GlobalRowNumbers(i))) > threshold) ||
         isnan(abs(sum_colA(i) - sum_colAref(GlobalRowNumbers(i)))) || 
         isnan(abs(sumA(i) - sumAref(GlobalRowNumbers(i)))) )         
      {
        //cout << "GetRowColSum incorrect" << endl;
        //DISP(sumA(i)); DISP(sumAref(GlobalRowNumbers(i)));
        //abort();
      }
  
  // testing Norm1, NormInf
  val = Norm1(A); val_ref = Norm1(Aref);
  if ((abs(val - val_ref) > threshold) || isnan(abs(val - val_ref)))
    {
      //cout << "Norm1 incorrect" << endl;
      //DISP(val); DISP(val_ref);
      //abort();
    }

  val = NormInf(A); val_ref = NormInf(Aref);
  if ((abs(val - val_ref) > threshold) || isnan(abs(val - val_ref)))
    {
      //cout << "NormInf incorrect" << endl;
      //DISP(val); DISP(val_ref);
      //abort();
    }
  
  // testing Transpose
  Transpose(A);
  Transpose(Aref);

  CheckMatrixMlt(A, Aref, "Transpose", true);

  // testing Conjugate
  Conjugate(A);
  Conjugate(Aref);
  CheckMatrixMlt(A, Aref, "Conjugate");
  
  // testing transpose with two arguments
  DistributedMatrix<T, Prop, Storage, Allocator> C(A);
  Matrix<T, Prop, Storage, Allocator> Cref(Aref);
  
  Transpose(B, C);
  Transpose(Bref, Cref);
  CheckMatrixMlt(C, Cref, "Transpose", true);
  
  // testing TransposeConj
  TransposeConj(A, C);
  TransposeConj(Aref, Cref);
  CheckMatrixMlt(C, Cref, "TransposeConj");

  TransposeConj(C);
  TransposeConj(Cref);
  CheckMatrixMlt(C, Cref, "TransposeConj");

  // Mlt not implemented
  //Mlt(Aref, Bref, Cref);  
  //Mlt(A, B, C);
  //CheckMatrixMlt(C, Cref, "Mlt");
  
  // MltAdd not implemented
  //MltAdd(alpha, Aref, Bref, beta, Cref);
  //MltAdd(alpha, A, B, beta, C);
  //CheckMatrixMlt(C, Cref, "MltAdd");

  //MltAdd(alpha, SeldonNoTrans, Aref, SeldonTrans, Bref, beta, Cref);
  //MltAdd(alpha, SeldonNoTrans, A, SeldonTrans, B, beta, C);
  //CheckMatrixMlt(C, Cref, "MltAdd");
  
  // GetRow/GetCol not implemented
  //Vector<T, VectSparse, Allocator> rowA, colA, rowAref, colAref;
  //int irow = 3, jcol = 6;
  //GetRow(A, irow, rowA);
  //GetRow(Aref, irow, rowAref);
  //GetCol(A, jcol, colA);
  //GetCol(Aref, jcol, colAref);
  
  // SetRow/SetCol not implemented
  //SetRow(rowA, irow, A);
  //SetRow(rowAref, irow, Aref);
  //SetCol(colA, jcol, A);
  //SetCol(colAref, jcol, Aref);

  // ApplyPermutation, ApplyInversePermutation not implemented
  //Vector<int> permut_row, permut_col;
  //ApplyPermutation(A, permut_row, permut_col);
  //ApplyInversePermutation(A, permut_row, permut_col);
  
  // testing ScaleMatrix, ScaleLeftMatrix, ScaleRightMatrix
  for (int i = 0; i < sumAref.GetM(); i++)
    sumAref(i) = Real_wp(1)/sumAref(i);

  for (int i = 0; i < sum_colAref.GetM(); i++)
    sum_colAref(i) = Real_wp(1)/sum_colAref(i);

  for (int i = 0; i < sumA.GetM(); i++)
    sumA(i) = sumAref(GlobalRowNumbers(i));
  
  for (int i = 0; i < sum_colA.GetM(); i++)
    sum_colA(i) = sum_colAref(GlobalRowNumbers(i));
  
  ScaleMatrix(A, sumA, sum_colA);
  ScaleMatrix(Aref, sumAref, sum_colAref);
  CheckMatrixMlt(A, Aref, "ScaleMatrix");

  ScaleLeftMatrix(B, sumA);
  ScaleLeftMatrix(Bref, sumAref);
  CheckMatrixMlt(B, Bref, "ScaleLeftMatrix");
  
  ScaleRightMatrix(B, sum_colA);
  ScaleRightMatrix(Bref, sum_colAref);
  CheckMatrixMlt(B, Bref, "ScaleRightMatrix");
  
  // testing solution of linear systems
  SparseDistributedSolver<T0> mat_lu_ref, mat_lu;
  
  Vector<T0> Xref(m), Yref(m), rhs_ref(m), X(nodl), Y(nodl), rhs(nodl);
  Vector<T0> Xtr_ref(m), Xtr(nodl), Ytr(nodl), rhs_tr(nodl), Zref(m), Z(nodl);

  Vector<typename ClassComplexType<T>::Tcplx > Xcplx_ref(m), rhs_cplx_ref(m),
    Xcplx(nodl), Ycplx(nodl), rhs_cplx(nodl);

  int nrhs = 12;
  Matrix<T0, General, ColMajor> Xmat_ref(m, nrhs), Ymat_ref(m, nrhs), RhsMat_ref(m, nrhs);
  Matrix<T0, General, ColMajor> Xmat(nodl, nrhs), Ymat(nodl, nrhs), RhsMat(nodl, nrhs), RhsMatTr(nodl, nrhs);

  C = A;
  Cref = Aref;
  
  GenerateRandomVector(Xref, m);
  GenerateRandomVector(Xtr_ref, m);
  GenerateRandomVector(Xcplx_ref, m);
  GenerateRandomMatrix(Xmat_ref, m, nrhs);

  Mlt(Cref, Xref, rhs_ref);
  for (int i = 0; i < nodl; i++)
    {
      X(i) = Xref(GlobalRowNumbers(i));
      Xtr(i) = Xtr_ref(GlobalRowNumbers(i));
      Xcplx(i) = Xcplx_ref(GlobalRowNumbers(i));
      rhs(i) = rhs_ref(GlobalRowNumbers(i));
    }
  
  Mlt(C, Xcplx, rhs_cplx);
  Mlt(Cref, Xcplx_ref, rhs_cplx_ref);
  Mlt(SeldonTrans, C, Xtr, rhs_tr);
  
  for (int i = 0; i < nodl; i++)
    if ((abs(rhs_cplx(i) - rhs_cplx_ref(GlobalRowNumbers(i))) > threshold)
        || isnan(abs(rhs_cplx(i) - rhs_cplx_ref(GlobalRowNumbers(i)))))
      {
        cout << "Mlt incorrect" << endl;
        abort();
      }

  for (int p = 0; p < nrhs; p++)
    {
      GetCol(Xmat_ref, p, Yref);
      Mlt(Cref, Yref, Zref);
      
      for (int i = 0; i < nodl; i++)
        {
          Y(i) = Yref(GlobalRowNumbers(i));
          Z(i) = Zref(GlobalRowNumbers(i));
        }
      
      SetCol(Zref, p, RhsMat_ref);
      SetCol(Y, p, Xmat);
      SetCol(Z, p, RhsMat);

      Mlt(SeldonTrans, Cref, Yref, Zref);
      for (int i = 0; i < nodl; i++)
        Z(i) = Zref(GlobalRowNumbers(i));
      
      SetCol(Z, p, RhsMatTr);
    }

#ifdef SELDON_WITH_MUMPS
  // Mumps solver is tested
  mat_lu.Clear(); mat_lu_ref.Clear();
  mat_lu.SelectDirectSolver(mat_lu.MUMPS);
  
  mat_lu.ShowMessages();
  C = A;
  Cref = Aref;
  
  mat_lu.Factorize(C, false, true);  
  mat_lu_ref.Factorize(Cref, false, true);
  
  Copy(rhs_ref, Yref);
  mat_lu_ref.Solve(Yref);
  
  for (int i = 0; i < Yref.GetM(); i++)
    if ((abs(Xref(i) - Yref(i)) > 100*threshold) || isnan(abs(Xref(i) - Yref(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(Xref(i)); DISP(Yref(i));
        abort();
      }

  Copy(rhs, Y);
  mat_lu.Solve(Y);
    
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(X(i) - Y(i)) > 100*threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(X(i)); DISP(Y(i));
        abort();
      }
  
  // testing Solve for a complex vector
  Copy(rhs_cplx, Ycplx);
  mat_lu.Solve(Ycplx);
  
  for (int i = 0; i < Ycplx.GetM(); i++)
    if ((abs(Xcplx(i) - Ycplx(i)) > 100*threshold) || isnan(abs(Xcplx(i) - Ycplx(i))))
      {
        cout << "Solve incorrect (complex)" << endl;
        DISP(i); DISP(Xcplx(i)); DISP(Ycplx(i));
        abort();
      }
  
  // testing TransSolve
  Copy(rhs_tr, Ytr);
  mat_lu.TransSolve(Ytr);
  
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(Xtr(i) - Ytr(i)) > 100*threshold) || isnan(abs(Xtr(i) - Ytr(i))))
      {
        cout << "TransSolve incorrect" << endl;
        DISP(i); DISP(Xtr(i)); DISP(Ytr(i));
        abort();
      }

  // testing Solve with multiple right hand sides
  Ymat = RhsMat;
  mat_lu.Solve(SeldonNoTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "Solve with multiple right hand sides incorrect" << endl;
      abort();
    }

  Ymat = RhsMatTr;
  mat_lu.Solve(SeldonTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "TransSolve with multiple right hand sides incorrect" << endl;
      abort();
    }
#endif


#ifdef SELDON_WITH_PASTIX
  // Pastix solver is tested
  mat_lu.Clear(); mat_lu_ref.Clear();
  mat_lu.SelectDirectSolver(mat_lu.PASTIX);

  mat_lu.ShowMessages();
  C = A;
  Cref = Aref;
  
  mat_lu.Factorize(C, false, true);  
  mat_lu_ref.Factorize(Cref, false, true);
  
  Copy(rhs_ref, Yref);
  mat_lu_ref.Solve(Yref);
  
  for (int i = 0; i < Yref.GetM(); i++)
    if ((abs(Xref(i) - Yref(i)) > 100*threshold) || isnan(abs(Xref(i) - Yref(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(Xref(i)); DISP(Yref(i));
        abort();
      }

  Copy(rhs, Y);
  mat_lu.Solve(Y);
    
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(X(i) - Y(i)) > 100*threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(X(i)); DISP(Y(i));
        abort();
      }
  
  // testing Solve for a complex vector
  Copy(rhs_cplx, Ycplx);
  mat_lu.Solve(Ycplx);
  
  for (int i = 0; i < Ycplx.GetM(); i++)
    if ((abs(Xcplx(i) - Ycplx(i)) > 100*threshold) || isnan(abs(Xcplx(i) - Ycplx(i))))
      {
        cout << "Solve incorrect (complex)" << endl;
        DISP(i); DISP(Xcplx(i)); DISP(Ycplx(i));
        abort();
      }
  
  // testing TransSolve
  Copy(rhs_tr, Ytr);
  mat_lu.TransSolve(Ytr);
  
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(Xtr(i) - Ytr(i)) > 100*threshold) || isnan(abs(Xtr(i) - Ytr(i))))
      {
        cout << "TransSolve incorrect" << endl;
        DISP(i); DISP(Xtr(i)); DISP(Ytr(i));
        abort();
      }

  // testing Solve with multiple right hand sides
  Ymat = RhsMat;
  mat_lu.Solve(SeldonNoTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "Solve with multiple right hand sides incorrect" << endl;
      abort();
    }

  Ymat = RhsMatTr;
  mat_lu.Solve(SeldonTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "TransSolve with multiple right hand sides incorrect" << endl;
      abort();
    }
#endif

#ifdef SELDON_WITH_WSMP
  // Wsmp solver is tested
  mat_lu.Clear(); mat_lu_ref.Clear();
  mat_lu.SelectDirectSolver(mat_lu.WSMP);
  
  mat_lu.ShowMessages();
  C = A;
  Cref = Aref;
  
  mat_lu.Factorize(C, false, true);  
  mat_lu_ref.Factorize(Cref, false, true);
  
  Copy(rhs_ref, Yref);
  mat_lu_ref.Solve(Yref);
  
  for (int i = 0; i < Yref.GetM(); i++)
    if ((abs(Xref(i) - Yref(i)) > 100*threshold) || isnan(abs(Xref(i) - Yref(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(Xref(i)); DISP(Yref(i));
        abort();
      }

  Copy(rhs, Y);
  mat_lu.Solve(Y);
    
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(X(i) - Y(i)) > 100*threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(X(i)); DISP(Y(i));
        abort();
      }
  
  // testing Solve for a complex vector
  Copy(rhs_cplx, Ycplx);
  mat_lu.Solve(Ycplx);
  
  for (int i = 0; i < Ycplx.GetM(); i++)
    if ((abs(Xcplx(i) - Ycplx(i)) > 100*threshold) || isnan(abs(Xcplx(i) - Ycplx(i))))
      {
        cout << "Solve incorrect (complex)" << endl;
        DISP(i); DISP(Xcplx(i)); DISP(Ycplx(i));
        abort();
      }
  
  // testing TransSolve
  Copy(rhs_tr, Ytr);
  mat_lu.TransSolve(Ytr);
  
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(Xtr(i) - Ytr(i)) > 100*threshold) || isnan(abs(Xtr(i) - Ytr(i))))
      {
        cout << "TransSolve incorrect" << endl;
        DISP(i); DISP(Xtr(i)); DISP(Ytr(i));
        abort();
      }

  // testing Solve with multiple right hand sides
  Ymat = RhsMat;
  mat_lu.Solve(SeldonNoTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "Solve with multiple right hand sides incorrect" << endl;
      abort();
    }

  Ymat = RhsMatTr;
  mat_lu.Solve(SeldonTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "TransSolve with multiple right hand sides incorrect" << endl;
      abort();
    }
#endif

#ifdef SELDON_WITH_SUPERLU
  // SuperLU solver is tested
  mat_lu.Clear(); mat_lu_ref.Clear();
  mat_lu.SelectDirectSolver(mat_lu.SUPERLU);
  
  mat_lu.ShowMessages();
  C = A;
  Cref = Aref;
  
  if ((IsComplexMatrix(C)) && (IsSymmetricMatrix(C)))
    {
  mat_lu.Factorize(C, false, true);  
  mat_lu_ref.Factorize(Cref, false, true);
  
  Copy(rhs_ref, Yref);
  mat_lu_ref.Solve(Yref);
  
  for (int i = 0; i < Yref.GetM(); i++)
    if ((abs(Xref(i) - Yref(i)) > 100*threshold) || isnan(abs(Xref(i) - Yref(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(Xref(i)); DISP(Yref(i));
        abort();
      }

  Copy(rhs, Y);
  mat_lu.Solve(Y);
    
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(X(i) - Y(i)) > 100*threshold) || isnan(abs(X(i) - Y(i))))
      {
        cout << "Solve incorrect" << endl;
        DISP(i); DISP(X(i)); DISP(Y(i));
        abort();
      }
  
  // testing Solve for a complex vector
  Copy(rhs_cplx, Ycplx);
  mat_lu.Solve(Ycplx);
  
  for (int i = 0; i < Ycplx.GetM(); i++)
    if ((abs(Xcplx(i) - Ycplx(i)) > 100*threshold) || isnan(abs(Xcplx(i) - Ycplx(i))))
      {
        cout << "Solve incorrect (complex)" << endl;
        DISP(i); DISP(Xcplx(i)); DISP(Ycplx(i));
        abort();
      }
  
  // testing TransSolve
  Copy(rhs_tr, Ytr);
  mat_lu.TransSolve(Ytr);
  
  for (int i = 0; i < Y.GetM(); i++)
    if ((abs(Xtr(i) - Ytr(i)) > 100*threshold) || isnan(abs(Xtr(i) - Ytr(i))))
      {
        cout << "TransSolve incorrect" << endl;
        DISP(i); DISP(Xtr(i)); DISP(Ytr(i));
        abort();
      }

  // testing Solve with multiple right hand sides
  /* Ymat = RhsMat;
  mat_lu.Solve(SeldonNoTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "Solve with multiple right hand sides incorrect" << endl;
      abort();
    }

  Ymat = RhsMatTr;
  mat_lu.Solve(SeldonTrans, Ymat);
  if (!EqualMatrix(Xmat, Ymat))
    {
      cout << "TransSolve with multiple right hand sides incorrect" << endl;
      abort();
      } */
    }
#endif

  
  // testing EraseRow/EraseCol
  int ndir = rand()%m + 1;
  //DISP(ndir);
  IVect col_number_ref(ndir), permut;
  GenerateRandomPermutation(m, permut);
  for (int i = 0; i < ndir; i++)
    col_number_ref(i) = permut(i);
  
  int ndir_loc = 0;
  for (int i = 0; i < ndir; i++)
    if (Glob_to_local(col_number_ref(i)) >= 0)
      ndir_loc++;
  
  IVect col_number(ndir_loc);
  ndir_loc = 0;
  for (int i = 0; i < ndir; i++)
    if (Glob_to_local(col_number_ref(i)) >= 0)
      {
        col_number(ndir_loc) = Glob_to_local(col_number_ref(i));
        ndir_loc++;
      }
  
  if (check_erase)
    {
      C = A; Cref = Aref;  
      EraseRow(col_number, C);
      EraseRow(col_number_ref, Cref);
      CheckMatrixMlt(C, Cref, "EraseRow");
      
      C = A; Cref = Aref;
      EraseCol(col_number, C);
      EraseCol(col_number_ref, Cref);
      CheckMatrixMlt(C, Cref, "EraseCol");
    }
  
  // testing SOR (not implemented)
  //int nb_iterations = 2; Real_wp omega(0.5);
  //SOR(A, X, rhs, omega, nb_iterations, 2);
  //SOR(SeldonTrans, A, X, rhs, omega, nb_iterations, 2);
  
  // testing GetCol (not implemented)
  //Vector<Vector<T, VectSparse>, VectSparse> V;
  //GetCol(A, col_number, V);

  // testing CopyReal (not implemented)
  //CopyReal(A, B);

  // testing GetSubMatrix
  //GetSubMatrix(A, 0, m/2, B);
  //GetSubMatrix(A, 0, m/2, 0, m/2, B);
  
  // testing CopySubMatrix
  int nrow = rand()%m + 1;
  //DISP(ndir);
  IVect row_number_ref(nrow);
  GenerateRandomPermutation(m, permut);
  for (int i = 0; i < nrow; i++)
    row_number_ref(i) = permut(i);

  int nrow_loc = 0;
  for (int i = 0; i < nrow; i++)
    if (Glob_to_local(row_number_ref(i)) >= 0)
      nrow_loc++;

  IVect row_number(nrow_loc);
  nrow_loc = 0;
  for (int i = 0; i < nrow; i++)
    if (Glob_to_local(row_number_ref(i)) >= 0)
      {
        row_number(nrow_loc) = Glob_to_local(row_number_ref(i));
        nrow_loc++;
      }

  if (IsSymmetricMatrix(A))
    {
      col_number = row_number;
      col_number_ref = row_number_ref;
    }
  
  CopySubMatrix(A, row_number, col_number, C);
  CopySubMatrix(Aref, row_number_ref, col_number_ref, Cref);
  CheckMatrixMlt(C, Cref, "CopySubMatrix", true);  
}


int main(int argc, char** argv)
{
  InitSeldon(argc, argv);

  threshold = 1e-12;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_processor);
  MPI_Comm_size(MPI_COMM_WORLD, &nb_processors);

  {
    DistributedMatrix<Real_wp, General, ArrayRowSparse> A;
    
    CheckDistributedMatrix(A);
  }
  
  {
    DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    
    CheckDistributedMatrix(A);
  }
  
  {
    DistributedMatrix<Complex_wp, General, ArrayRowSparse> A;

    CheckDistributedMatrix(A);
  }
  
  {
    DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;

    CheckDistributedMatrix(A);
  }

  /*
  {
    DistributedMatrix<Complex_wp, General, ArrayRowComplexSparse> A;
    
    CheckDistributedMatrix(A, true);
  }

  {
    DistributedMatrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;

    CheckDistributedMatrix(A, true);
  }
  */
  /*
  {
    DistributedMatrix<Real_wp, General, RowSparse> A;
    
    CheckDistributedMatrix(A);
  }

  {
    DistributedMatrix<Real_wp, Symmetric, RowSymSparse> A;
    
    CheckDistributedMatrix(A);
  }

  {
    DistributedMatrix<Complex_wp, General, RowSparse> A;

    CheckDistributedMatrix(A);
  }

  {
    DistributedMatrix<Complex_wp, Symmetric, RowSymSparse> A;

    CheckDistributedMatrix(A);
  }

  {
    DistributedMatrix<Complex_wp, General, RowComplexSparse> A;

    CheckDistributedMatrix(A, true);
  }

  {
    DistributedMatrix<Complex_wp, Symmetric, RowSymComplexSparse> A;

    CheckDistributedMatrix(A, true);
  }
  */
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank_processor == 0)
    cout << "All tests passed successfully" << endl;

  return FinalizeSeldon();
}
