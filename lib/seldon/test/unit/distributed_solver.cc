#include "SeldonLib.hxx"

using namespace Seldon;

int rank_processor, nb_processors, root_processor;

#include "Poisson.cxx"

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);

  {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_processor);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_processors);
    root_processor = 0;
  
  double w = 0.31*M_PI, w2 = w*w;
  double L = 10.0;
  int nb_points_glob = 100;
  int Nglob = nb_points_glob*nb_points_glob;
  double dx = L/(nb_points_glob + 1);
  


  // splitting the 2-D domain between processors
  int nb_proc_row = 2, nb_proc_col = 2;
  if (nb_proc_row*nb_proc_col != nb_processors)
    {
      cout << "Run it with the correct number of processors" << endl;
      abort();
    }

  IVect GlobalRowNumbers;
  IVect ListProcInterac;
  Vector<IVect> ListRowInterac;
  IVect OverlapRowNumbers, OriginalProcNumbers;
  MPI_Comm comm = MPI_COMM_WORLD;

  // function to distribute dofs between processors
  // code in Poisson.cxx
  int i0_proc, i1_proc, j0_proc, j1_proc;
  DistributePoissonDofs(nb_proc_row, nb_proc_col, nb_points_glob,
			comm, GlobalRowNumbers, ListProcInterac,
			ListRowInterac, OverlapRowNumbers, OriginalProcNumbers,
			i0_proc, i1_proc, j0_proc, j1_proc);

  int N = GlobalRowNumbers.GetM();
  int Nx = i1_proc-i0_proc; // int Ny = j1_proc-j0_proc;
  
  // computes Helmholtz matrix
  DistributedMatrix<double, General, ArrayRowSparse> A;
  A.Reallocate(N, N);
  A.Init(Nglob, &GlobalRowNumbers, &OverlapRowNumbers, &OriginalProcNumbers,
	 N, 1, &ListProcInterac, &ListRowInterac, comm);
  
  ComputeHelmholtzMatrix(A, dx, i1_proc-i0_proc, j1_proc-j0_proc, w2);
  A.WriteText("mat.dat");
  
  // factorizes matrix
  bool iterative = true;
#ifdef SELDON_WITH_MUMPS
  iterative = false;
#endif
  SparseDistributedSolver<double> mat_lu;
  
  if (!iterative)
    mat_lu.Factorize(A);
  
  // computes right hand side (dirac at the center)
  Vector<double> x_sol(N);
  x_sol.Fill(0.0);
  int iglob0 = nb_points_glob/2;
  if ((iglob0 >= i0_proc) && (iglob0 < i1_proc))
    if ((iglob0 >= j0_proc) && (iglob0 < j1_proc))
      {
	int i = iglob0 - i0_proc;
	int j = iglob0 - j0_proc;
	x_sol(j*Nx + i) = 1.0;
      }
  
  // solves linear sytem  
  if (!iterative)
    mat_lu.Solve(x_sol);
  else
    {
      DistributedVector<double> source(OverlapRowNumbers, comm), sol(OverlapRowNumbers, comm);
      
      source.SetData(x_sol);
      sol.Reallocate(N);

      Preconditioner_Base<double> prec;
      Iteration<double> iter(10000, 1e-6);
      iter.HideMessages();
      if (rank_processor == 0)
        iter.ShowFullHistory();
      
      Cg(A, sol, source, prec, iter);
      for (int i = 0; i < N; i++)
        x_sol(i) = sol(i);

      source.Nullify();
    }
  
  // writes solution
  if (rank_processor != 0)
    {
      IVect num(4);
      num(0) = i0_proc; num(1) = i1_proc; num(2) = j0_proc; num(3) = j1_proc;
      MPI_Send(num.GetData(), 4, MPI_INTEGER, 0, 2, comm);
      
      MPI_Send(x_sol.GetData(), x_sol.GetM(), MPI_DOUBLE, 0, 3, comm);
    }
  else
    {
      Matrix<double> x_mat(nb_points_glob, nb_points_glob);
      for (int i = i0_proc; i < i1_proc; i++)
	for (int j = j0_proc; j < j1_proc; j++)
	  x_mat(i, j) = x_sol((j-j0_proc)*Nx + i-i0_proc);
      
      IVect num(4);
      Vector<double> xtmp; MPI_Status status;
      for (int rank = 1; rank < nb_processors; rank++)
	{
	  MPI_Recv(num.GetData(), 4, MPI_INTEGER, rank, 2, comm, &status);
	  int i0 = num(0), i1 = num(1), j0 = num(2), j1 = num(3);
	  xtmp.Reallocate((i1-i0)*(j1-j0));
	  MPI_Recv(xtmp.GetData(), xtmp.GetM(), MPI_DOUBLE, rank, 3, comm, &status);
	  
	  for (int i = i0; i < i1; i++)
	    for (int j = j0; j < j1; j++)
	      x_mat(i, j) = xtmp((j-j0)*(i1-i0) + i-i0);
	}
      
  
      x_mat.WriteText("Sol.dat");
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "Destroying variables" << endl;
    } 
  
  return FinalizeSeldon();
}
