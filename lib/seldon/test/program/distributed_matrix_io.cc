#include "SeldonLib.hxx"

using namespace Seldon;

int main(int argc, char** argv)
{
  InitSeldon(argc, argv);

  if (argc != 3)
    {
      cout << "Please provide the file where the matrix and right handside is stored " << endl;
      cout << "Usage : ./exe matrix.dat rhs.dat" << endl;
      abort();
    }

  {
    typedef complex<double> Complexe;
    
  string file_matrix(argv[1]);
  string file_rhs(argv[2]);
  
  MPI_Comm comm = MPI_COMM_WORLD;

  int rank_proc; MPI_Comm_rank(comm, &rank_proc);
  int nb_proc; MPI_Comm_size(comm, &nb_proc);

  // global matrix A is read in processor 0
  // num_proc(i) contains the processor number for row i
  Matrix<Complexe, Symmetric, ArrayRowSymSparse> A;
  Vector<int> num_proc; int nb_rows = 0;
  if (rank_proc == 0)
    {
      A.ReadText(file_matrix);
      nb_rows = A.GetM();
      num_proc.Reallocate(nb_rows);
      num_proc.Zero();

      // simple distribution of rows
      // processor 0 receives rows 0, 1, .. N0
      // processor 1 receives rows N0+1, N0+2, ... N1
      // processor 2 receives rows N1+1, ... N2
      if (nb_proc > 1)
	{
	  int nb_rows_div = int(ceil(double(nb_rows)/nb_proc));
	  int p = nb_rows - (nb_rows_div-1)*nb_proc; 
	  
	  int offset = 0;
	  for (int i = 0; i < nb_proc; i++)
	    {
	      int nb_rows_p = nb_rows_div-1;
	      if (i < p)
		nb_rows_p = nb_rows_div;
		
	      for (int j = 0; j < nb_rows_p; j++)
		num_proc(offset+j) = i;
	      
	      offset += nb_rows_p;
	    }
	}

      cout << "Chosen distribution of rows = " << num_proc << endl;
    }
  
  // the matrix is distributed between processors
  DistributedMatrix<Complexe, Symmetric, ArrayRowSymSparse> B;
  DistributedMatrixIntegerArray Binfo;
  DistributeCentralizedMatrix(A, comm, 0, num_proc, Binfo, B);
  
  B.WriteText("mat_dist.dat");

  // then the matrix is factorized
  SparseDistributedSolver<Complexe> mat_lu;

  int m = B.GetM();
  mat_lu.Factorize(B);

  // the right hand side is read and solved
  Vector<Complexe> x, xloc(m);
  x.Read(file_rhs);
  for (int i = 0; i < Binfo.GlobalRowNumbers.GetM(); i++)
    xloc(i) = x(Binfo.GlobalRowNumbers(i));  
  
  mat_lu.Solve(xloc);

  // then we centralize the solution (not optimal here)
  // for a gathering more optimal see Mumps.cxx
  x.Zero();
  for (int i = 0; i < Binfo.GlobalRowNumbers.GetM(); i++)
    x(Binfo.GlobalRowNumbers(i)) = xloc(i);  
  
  Vector<Complexe> y(x.GetM()); Vector<int64_t> xtmp;
  MpiReduce(comm, x, xtmp, y, x.GetM(), MPI_SUM, 0);
  
  if (rank_proc == 0)
    y.Write("sol.dat");

  }
  
  return FinalizeSeldon();
}
