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
  
  double w = 0.01*M_PI, w2 = w*w;
  double L = 10.0;
  int nb_points_glob = 100;
  int Nglob = nb_points_glob*nb_points_glob;
  double dx = L/(nb_points_glob - 1);

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
  
  // computes Laplacian matrix
  DistributedMatrix<double, Symmetric, ArrayRowSymSparse> A;
  A.Reallocate(N, N);
  A.Init(Nglob, &GlobalRowNumbers, &OverlapRowNumbers, &OriginalProcNumbers,
	 N, 1, &ListProcInterac, &ListRowInterac, comm);
  
  ComputeHelmholtzMatrix(A, dx, i1_proc-i0_proc, j1_proc-j0_proc, 0.0);
  A.WriteText("mat.dat");

  // object for computing eigenvalues
  SparseEigenProblem<double, DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse>,
		     DistributedMatrix<Real_wp, Symmetric, ArrayRowSymSparse> > var_eig;

  // int type_solver = TypeEigenvalueSolver::SLEPC;
  
  var_eig.SetStoppingCriterion(1e-12);
  var_eig.SetNbAskedEigenvalues(10);
  var_eig.SetComputationalMode(var_eig.SHIFTED_MODE);
  var_eig.SetNbMaximumIterations(20000);
  var_eig.SetPrintLevel(1);

  var_eig.InitMatrix(A);

  var_eig.SetTypeSpectrum(var_eig.CENTERED_EIGENVALUES, w2, var_eig.SORTED_MODULUS);

  Vector<double> lambda, lambda_imag;
  Matrix<double, General, ColMajor> eigen_vec;

  // manual selection of the solver
  //GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec, type_solver);

  // automatic selection of the solver
  GetEigenvaluesEigenvectors(var_eig, lambda, lambda_imag, eigen_vec);

  if (rank_processor == 0)
    DISP(lambda);
  
  Vector<double> x_sol(N);
  for (int p = 0; p < eigen_vec.GetN(); p++)
    {
      GetCol(eigen_vec, p, x_sol);
      
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
	  
	  
	  x_mat.WriteText(string("Mode" + to_str(p) + ".dat"));
	}
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  cout << "Destroying variables" << endl;
    } 
  
  return FinalizeSeldon();
}
