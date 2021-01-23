#include "SeldonLib.hxx"

using namespace Seldon;

Real_wp threshold;

template<class T>
void GetRand(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T>& x)
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
    GetRand(x(i));
}

template<class T>
bool CheckVector(Vector<T>& x)
{
  bool test = true;
  T icplx, coef;
  SetComplexReal(x.GetM(), icplx);
  SetComplexOne(coef); coef /= icplx;
  for (int i = 0; i < x.GetM(); i++)
    {
      SetComplexReal(i, icplx);
      if ((abs(x(i) - icplx)*coef > threshold) || isnan(x(i)))
        test = false;
    }

  return test;
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
    if ((abs(x(i) - y(i)) > eps) || isnan(x(i)) || isnan(y(i)))
      return false;
  
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
  sum = zero; Vector<Real_wp> coef(n);
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

void FindPartitioningDof(int m, Vector<IVect>& list_proc,
			 IVect& Glob_to_local, IVect& GlobalRowNumbers,
			 IVect& OverlapRowNumbers, IVect& OverlapProcNumbers,
			 int& nodl_scalar, IVect& ProcSharingRows,
			 Vector<IVect>& SharingRowNumbers,
			 Vector<IVect>& all_glob, MPI_Comm comm)
{  
  int nb_processors; MPI_Comm_size(comm, &nb_processors);
  int rank_processor; MPI_Comm_rank(comm, &rank_processor);
  
  // for each dof determining processors that share this dof
  list_proc.Reallocate(m);
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
  Glob_to_local.Reallocate(m); Glob_to_local.Fill(-1);
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
  
  GlobalRowNumbers.Reallocate(nodl);
  for (int i = 0; i < m; i++)
    if (Glob_to_local(i) >= 0)
      GlobalRowNumbers(Glob_to_local(i)) = i;
  
  //DISP(nodl); DISP(GlobalRowNumbers);
  
  // then OverlapRowNumbers/OverlapProcNumbers
  OverlapRowNumbers.Reallocate(nodl_overlap);
  OverlapProcNumbers.Reallocate(nodl_overlap);
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
  nodl_scalar = nodl; // nb_u = 1
  int nb_proc = 0;
  IVect IndexProc(nb_processors); IndexProc.Fill(-1);
  for (int p = 0; p < nb_processors; p++)
    if (NbSharedDofPerProc(p) > 0)
      {
        IndexProc(p) = nb_proc;
        nb_proc++;
      }
  
  ProcSharingRows.Reallocate(nb_proc);
  SharingRowNumbers.Reallocate(nb_proc);
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

  if (rank_processor == 0)
    {
      all_glob.Reallocate(nb_processors);
      all_glob(0) = GlobalRowNumbers;
      MPI_Status status;
      for (int p = 1; p < nb_processors; p++)
	{
	  int nb_dof_p = 0;
	  MPI_Recv(&nb_dof_p, 1, MPI_INTEGER, p, 23, comm, &status);
	  all_glob(p).Reallocate(nb_dof_p);
	  MPI_Recv(all_glob(p).GetData(), nb_dof_p, MPI_INTEGER, p, 24, comm, &status);

	  //DISP(p); DISP(all_glob);
	}
    }
  else
    {
      MPI_Send(&nodl, 1, MPI_INTEGER, 0, 23, comm);
      MPI_Send(GlobalRowNumbers.GetData(), nodl, MPI_INTEGER, 0, 24, comm);
    }
  
}


template<class MatrixSeq, class MatrixPar>
void DistributeMatrixProcessor(const MatrixSeq& Aref, MatrixPar& A, MPI_Comm comm,
                               const Vector<IVect>& list_proc, const IVect& Glob_to_local)
{
  int rank_processor; MPI_Comm_rank(comm, &rank_processor);
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

template<class T, class Prop, class Storage, class Allocator>
void CheckSparseCholesky(Matrix<T, Prop, Storage, Allocator>& Aref, MPI_Comm comm)
{
  
  // symmetric matrix is read in a file
  Aref.ReadText("test/matrix/MhSparse.dat");

  // rows are distributed among processors
  int m = Aref.GetM(); Vector<IVect> list_proc;
  IVect Glob_to_local, GlobalRowNumbers;
  IVect OverlapRowNumbers, OverlapProcNumbers;
  int nodl_scalar; IVect ProcSharingRows;
  Vector<IVect> SharingRowNumbers, all_glob;
  FindPartitioningDof(m, list_proc, Glob_to_local, GlobalRowNumbers, OverlapRowNumbers,
		      OverlapProcNumbers, nodl_scalar, ProcSharingRows, SharingRowNumbers,
		      all_glob, comm);
  
  // distributing Aref to all the processors
  DistributedMatrix<T, Prop, Storage, Allocator> A;
  int nodl = nodl_scalar;
  A.Reallocate(nodl, nodl);
  A.Init(m, &GlobalRowNumbers, &OverlapRowNumbers, &OverlapProcNumbers,
         nodl_scalar, 1, &ProcSharingRows, &SharingRowNumbers, comm);

  DistributeMatrixProcessor(Aref, A, comm, list_proc, Glob_to_local);

  // creation of a right hand side b = A*[0;1;...;n-1]
  Vector<T> xref(Aref.GetM()), bref(Aref.GetM()), yref(Aref.GetM());
  xref.Fill();
  Mlt(Aref, xref, bref);
  
  /*********************
   * ArrayRowSymSparse *
   *********************/
  
  // testing SparseCholeskySolver
  SparseCholeskySolver<T> mat_lu_ref;
  DistributedCholeskySolver<T> mat_lu;
  
  mat_lu_ref.SelectDirectSolver(mat_lu_ref.PASTIX);
  mat_lu_ref.Factorize(Aref);

  yref.Fill();
  xref = yref;
  mat_lu_ref.Mlt(SeldonTrans, xref);
  mat_lu_ref.Mlt(SeldonNoTrans, xref);

  mat_lu_ref.Solve(SeldonNoTrans, xref);
  mat_lu_ref.Solve(SeldonTrans, xref);
  
  if (!CheckVector(xref))
    {
      cout << "MltCholesky incorrect" << endl;
      abort();
    }

  xref = bref;
  mat_lu_ref.Solve(SeldonNoTrans, xref);
  mat_lu_ref.Solve(SeldonTrans, xref);

  if (!CheckVector(xref))
    {
      cout << "SolveCholesky incorrect" << endl;
      abort();
    }

  int nb_processors; MPI_Comm_size(comm, &nb_processors);
  int rank_processor; MPI_Comm_rank(comm, &rank_processor);
  if (nb_processors >= 1)
    {
      Vector<T> x(nodl), y(nodl), b(nodl), xt(nodl);
      for (int i = 0; i < nodl; i++)
	{
	  b(i) = bref(GlobalRowNumbers(i));
	  xt(i) = GlobalRowNumbers(i);
	}
      
      mat_lu.SelectDirectSolver(mat_lu.PASTIX);
      mat_lu.Factorize(A);

      x = xt;
      mat_lu.Mlt(SeldonNoTrans, x);
      mat_lu.Solve(SeldonNoTrans, x);
            
      if (!EqualVector(x, xt))
	{
	  cout << "MltCholesky incorrect" << endl;
	  abort();
	}

      x = xt;
      mat_lu.Mlt(SeldonTrans, x);
      mat_lu.Solve(SeldonTrans, x);

      //DISP(x);
      if (!EqualVector(x, xt))
	{
	  cout << "MltCholesky incorrect" << endl;
	  abort();
	}

      x = b;
      mat_lu.Solve(SeldonNoTrans, x);
      mat_lu.Solve(SeldonTrans, x);

      MPI_Barrier(comm);
      if (!EqualVector(x, xt))
	{
	  cout << "SolveCholesky incorrect" << endl;
	  abort();
	}

      Matrix<T, General, ArrayRowSparse> L(m, m);
      for (int i = 0; i < m; i++)
	{
	  x.Zero();
	  for (int j = 0; j < nodl; j++)
	    if (GlobalRowNumbers(j) == i)
	      x(j) = 1.0;

	  mat_lu.Mlt(SeldonNoTrans, x);

	  if (rank_processor == 0)
	    {
	      MPI_Status status;
	      for (int p = 0; p < nb_processors; p++)
		{
		  Vector<T> xp(all_glob(p).GetM());
		  if (p == 0)
		    xp = x;
		  else
		    MPI_Recv(xp.GetData(), xp.GetM(), MPI_DOUBLE, p, 25, comm, &status);

		  for (int j = 0; j < xp.GetM(); j++)
		    if (abs(xp(j)) != T(0))
		      L.AddInteraction(all_glob(p)(j), i, xp(j));
		}
	    }
	  else
	    {
	      MPI_Send(x.GetData(), nodl, MPI_DOUBLE, 0, 25, comm);
	    }
	}

      if (rank_processor == 0)
	L.WriteText("L.dat");

      Matrix<T, General, ArrayRowSparse> Lt(m, m);
      for (int i = 0; i < m; i++)
	{
	  x.Zero();
	  for (int j = 0; j < nodl; j++)
	    if (GlobalRowNumbers(j) == i)
	      x(j) = 1.0;

	  //DISP(i); DISP(x);
	  mat_lu.Mlt(SeldonTrans, x);
	  //DISP(x);
	  
	  if (rank_processor == 0)
	    {
	      MPI_Status status;
	      for (int p = 0; p < nb_processors; p++)
		{
		  Vector<T> xp(all_glob(p).GetM());
		  if (p == 0)
		    xp = x;
		  else
		    MPI_Recv(xp.GetData(), xp.GetM(), MPI_DOUBLE, p, 25, comm, &status);

		  for (int j = 0; j < xp.GetM(); j++)
		    if (abs(xp(j)) != T(0))
		      Lt.AddInteraction(all_glob(p)(j), i, xp(j));
		}
	    }
	  else
	    {
	      MPI_Send(x.GetData(), nodl, MPI_DOUBLE, 0, 25, comm);
	    }
	}

      if (rank_processor == 0)
	Lt.WriteText("Lt.dat");

      MPI_Barrier(comm);
      if (rank_processor == 0)
	{
	  Transpose(Lt);
	  Add(T(-1), Lt, L); DISP(MaxAbs(L));
	  if ((isnan(MaxAbs(L))) || (MaxAbs(L) > threshold))
	    {
	      cout << "Mlt incorrect" << endl;
	      abort();
	    }
	}
    }
}

void SplitWorldCommunicator(MPI_Comm& comm, string& suffix)
{
  if (true)
    {
      // no split
      comm = MPI_COMM_WORLD;
      suffix = string("");
    }
  else
    {
      // splitting in two communicators
      int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      int nb_proc; MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
      if (nb_proc < 4)
	{
	  cout << "This unitary test must be launched on four processors at least" << endl;
	  abort();
	}
      
      int color = 0;
      int key = rank;
      suffix = "_C0";
      if (rank >= nb_proc/2)
	{
	  color = 1;
	  key -= nb_proc/2;
	  suffix = "_C1";
	}
      
      MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm);
    }  
}

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);
  
  //srand(time(NULL));
  threshold = 1e-12;

  MPI_Comm comm; string suffix;
  SplitWorldCommunicator(comm, suffix);

  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSparseCholesky(A, comm);
  }

  MPI_Barrier(comm);

  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
    cout << "All tests passed successfully" << endl;
    
  return FinalizeSeldon();
}
