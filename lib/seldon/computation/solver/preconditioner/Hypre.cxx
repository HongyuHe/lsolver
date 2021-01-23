#ifndef SELDON_FILE_HYPRE_CXX

#include "Hypre.hxx"
#include "HypreInline.cxx"

namespace Seldon
{

  //! default constructor
  template<class T>
  HyprePreconditioner<T>::HyprePreconditioner()
  {
    solver_created = false;
    nodl_scalar = 0;
    nb_u = 1;

    comm = MPI_COMM_SELF;    
    print_level = 0;

    type_preconditioner = BOOMER_AMG;

    ProcNumber = NULL;
    DofNumber = NULL;
    
    amg_max_levels = 20;
    amg_num_sweeps = 1;
    amg_smoother = HYBRID_GS_BACKWARD;

    sai_filter = 0.1;
    sai_threshold = 0.1;
    sai_max_levels = 1;
    sai_sym = 1;

    euclid_use_ilut = false;
    euclid_level = 0;
    euclid_threshold = 0;
    euclid_droptol = 0.01;
  }


  //! destructor
  template<class T>
  HyprePreconditioner<T>::~HyprePreconditioner()
  {
    Clear();
  }


  //! erases the current object
  template<class T>
  void HyprePreconditioner<T>::Clear()
  {
    if (solver_created)
      {
	if (type_preconditioner == BOOMER_AMG)
	  HYPRE_BoomerAMGDestroy(solver);
	else if (type_preconditioner == PARASAILS)
	  HYPRE_ParaSailsDestroy(solver);
	else if (type_preconditioner == EUCLID)
	  HYPRE_EuclidDestroy(solver);
	else if (type_preconditioner == AMS)
	  HYPRE_AMSDestroy(solver);

	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(vec_b);
	HYPRE_IJVectorDestroy(vec_x);
      }

    solver_created = false;
  }


  //! sets parameters from a keyword and associated parameters
  template<class T>
  void HyprePreconditioner<T>::SetInputPreconditioning(const string& keyword, const Vector<string>& params)
  {
    if (keyword == "AMG")
      {
	type_preconditioner = BOOMER_AMG;
	for (int k = 0; k < params.GetM(); k += 2)
	  {
	    if (params(k) == "MaximumLevel")
	      amg_max_levels = to_num<int>(params(k+1));
	    else if (params(k) == "NumberSweeps")
	      amg_num_sweeps = to_num<int>(params(k+1));
	    else if (params(k) == "Smoother")
	      {
		if (params(k+1) == "Jacobi")
		  amg_smoother = JACOBI;
		else if (params(k+1) == "GaussSeidel")
		  amg_smoother = GS_PAR_SEQ;
		else if (params(k+1) == "HybridGaussSeidel")
		  amg_smoother = HYBRID_GS_BACKWARD;
		else if (params(k+1) == "HybridGaussSeidelForward")
		  amg_smoother = HYBRID_GS_FORWARD;
		else if (params(k+1) == "HybridGaussSeidelSymmetric")
		  amg_smoother = HYBRID_GS_SYMMETRIC;
		else if (params(k+1) == "L1_GaussSeidel")
		  amg_smoother = L1_GAUSS_SEIDEL;
		else if (params(k+1) == "Chebyshev")
		  amg_smoother = CHEBYSHEV;
		else if (params(k+1) == "FCF_Jacobi")
		  amg_smoother = FCF_JACOBI;
		else if (params(k+1) == "L1_Jacobi")
		  amg_smoother = L1_JACOBI;
		else
		  amg_smoother = to_num<int>(params(k+1));
	      }
	  }
      }
    else if (keyword == "ParaSails")
      {
	type_preconditioner = PARASAILS;
	for (int k = 0; k < params.GetM(); k += 2)
	  {
	    if (params(k) == "Filter")
	      sai_filter = to_num<HYPRE_Real>(params(k+1));
	    else if (params(k) == "Threshold")
	      sai_threshold = to_num<HYPRE_Real>(params(k+1));
	    else if (params(k) == "MaximumLevel")
	      sai_max_levels = to_num<int>(params(k+1));
	    else if (params(k) == "Symmetry")
	      sai_sym = to_num<int>(params(k+1));
	  }
      }
    else if (keyword == "Euclid")
      {
	type_preconditioner = EUCLID;
	for (int k = 0; k < params.GetM(); k += 2)
	  {
	    if (params(k) == "Threshold")
	      euclid_threshold = to_num<HYPRE_Real>(params(k+1));
	    else if (params(k) == "Droptol")
	      euclid_droptol = to_num<HYPRE_Real>(params(k+1));
	    else if (params(k) == "Level")
	      euclid_level = to_num<int>(params(k+1));
	    else if (params(k) == "Algorithm")
	      {
		if (params(k+1) == "ILU")
		  euclid_use_ilut = false;
		else if (params(k+1) == "ILUT")
		  euclid_use_ilut = true;
	      }
	  }
      }
    else if (keyword == "AMS")
      {
	type_preconditioner = AMS;
      }
    else
      {
	cout << "Unknown preconditioning : " << keyword << endl;
	abort();
      }
  }

  
  //! main function constructing hypre preconditioning from a given sparse matrix
  template<class T> template<class Prop, class Storage, class Allocator>
  void HyprePreconditioner<T>
  ::ConstructPreconditioner(DistributedMatrix<HYPRE_Complex, Prop, Storage, Allocator>& A0,
			    bool keep_matrix)
  {
    // previous preconditioning is cleared if existing
    Clear();
	
    comm = A0.GetCommunicator();
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    General prop; Vector<int> row_numbers;
    Vector<int> size_rows, cols; Vector<HYPRE_Complex> values; 

    if (nb_proc > 1)
      {
	ProcNumber = &A0.GetProcessorSharingRows();
	DofNumber = &A0.GetSharingRowNumbers();
	nodl_scalar = A0.GetNodlScalar();
	nb_u = A0.GetNbScalarUnknowns();
        
	if (keep_matrix)
	  {
	    DistributedMatrix<HYPRE_Complex, Prop, Storage, Allocator> A(A0);
	    AssembleDistributed(A, prop, comm, row_numbers,
				local_rows, size_rows, cols, values, false, true);
	  }
	else
	  AssembleDistributed(A0, prop, comm, row_numbers,
			      local_rows, size_rows, cols, values, false, true);
      }
    else
      {
	int N = A0.GetM();
	ConvertToCSR(A0, prop, size_rows, cols, values);
	if (!keep_matrix)
	  A0.Clear();

	row_numbers.Reallocate(N);
	row_numbers.Fill();
	local_rows.Reallocate(N);
	local_rows.Fill();
      }
    
    // we call a non-template function to avoid many instantiations of a lengthy code
    FinalizePreconditioner(row_numbers, size_rows, cols, values);
  }

  
  //! internal function to finalize the computation of preconditioning
  template<class T>
  void HyprePreconditioner<T>
  ::FinalizePreconditioner(Vector<int>& row_numbers, Vector<int>& size_rows,
			   Vector<int>& cols, Vector<HYPRE_Complex>& values)
  {    
    int Nloc = local_rows.GetM();
    int ilower = row_numbers(0);
    int iupper = row_numbers(Nloc-1);
      
    // size_rows is a ptr array (incremented numbers)
    // we decrement to have the size of rows
    for (int i = 0; i < Nloc; i++)
      size_rows(i) = size_rows(i+1) - size_rows(i);

    // Create the matrix.
    // Note that this is a square matrix, so we indicate the row partition
    // size twice (since number of rows = number of cols)
    HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);
  
    // Choose a parallel csr format storage
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    
    // Initialize before setting coefficients
    HYPRE_IJMatrixInitialize(A);
    
    // setting coefficients
    HYPRE_IJMatrixSetValues(A, Nloc, size_rows.GetData(), row_numbers.GetData(), cols.GetData(), values.GetData());

    // we clear intermediate arrays cols, values and size_rows
    cols.Clear(); values.Clear(); size_rows.Clear();

    // Assemble after setting the coefficients
    HYPRE_IJMatrixAssemble(A);
    
    // Get the parcsr matrix object to use
    HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
    
    // create two vectors for x and b
    HYPRE_IJVectorCreate(comm, ilower, iupper, &vec_b);
    HYPRE_IJVectorSetObjectType(vec_b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(vec_b);
    
    HYPRE_IJVectorCreate(comm, ilower, iupper, &vec_x);
    HYPRE_IJVectorSetObjectType(vec_x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(vec_x);

    HYPRE_IJVectorGetObject(vec_b, (void **) &par_b);
    HYPRE_IJVectorGetObject(vec_x, (void **) &par_x);

    // now constructing the preconditioning
    solver_created = true;
    if (type_preconditioner == BOOMER_AMG)
      {
	HYPRE_BoomerAMGCreate(&solver);	
	HYPRE_BoomerAMGSetPrintLevel(solver, print_level);
	
	HYPRE_BoomerAMGSetOldDefault(solver); // Falgout coarsening with modified classical interpolaiton
	HYPRE_BoomerAMGSetRelaxType(solver, amg_smoother);
	HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   // uses C/F relaxation
	HYPRE_BoomerAMGSetNumSweeps(solver, amg_num_sweeps);   // Sweeps on each level
	HYPRE_BoomerAMGSetMaxLevels(solver, amg_max_levels);  // maximum number of levels 
	HYPRE_BoomerAMGSetTol(solver, 0);      // conv. tolerance
	HYPRE_BoomerAMGSetMaxIter(solver, 1); // only one iteration for preconditioning
	
	// Now setup
	HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);

	// to avoid printing each time Solve is called
	HYPRE_BoomerAMGSetPrintLevel(solver, 0);
      }
    else if (type_preconditioner == PARASAILS)
      {
	// Now set up the ParaSails preconditioner and specify any parameters
	HYPRE_ParaSailsCreate(comm, &solver);
	
	// Set some parameters (See Reference Manual for more parameters)
	HYPRE_ParaSailsSetParams(solver, sai_threshold, sai_max_levels);
	HYPRE_ParaSailsSetFilter(solver, sai_filter);
	HYPRE_ParaSailsSetSym(solver, sai_sym);
	HYPRE_ParaSailsSetLogging(solver, print_level);

	// setup
	HYPRE_ParaSailsSetup(solver, parcsr_A, par_b, par_x);
      }
    else if (type_preconditioner == EUCLID)
      {
	HYPRE_EuclidCreate(comm, &solver);
	HYPRE_EuclidSetLevel(solver, euclid_level);

	HYPRE_EuclidSetSparseA(solver, euclid_threshold);
	if (euclid_use_ilut)
	  HYPRE_EuclidSetILUT(solver, euclid_droptol);
	
	if (print_level > 0)
	  {
	    HYPRE_EuclidSetStats(solver, 1);
	    HYPRE_EuclidSetMem(solver, 1);
	  }

	// setup
	HYPRE_EuclidSetup(solver, parcsr_A, par_b, par_x);
      }
    else
      {
	cout << "Unknown preconditioning : " << type_preconditioner << endl;
	abort();
      }
  }
  

  //! applies preconditioning z = M r (or its transpose z = M^T r)
  template<class T>
  void HyprePreconditioner<T>
  ::Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	  const Vector<T>& r, Vector<T>& z)
  {
    HYPRE_Complex* b_data = par_b->local_vector->data;
    HYPRE_Complex* x_data = par_x->local_vector->data;
    Vector<int>& num = local_rows;

    // we extract the original values of r to put in par_b (shared values are ignored)
    for (int i = 0; i < num.GetM(); i++)
      b_data[i] = r(num(i));

    for (int i = 0; i < num.GetM(); i++)
      x_data[i] = 0;

    // the appropriate hypre solver is called
    if (trans.NoTrans() || A.IsSymmetric())
      {
	switch(type_preconditioner)
	  {
	  case BOOMER_AMG:
	    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);
	    break;
	  case PARASAILS:
	    HYPRE_ParaSailsSolve(solver, parcsr_A, par_b, par_x);
	    break;
	  case EUCLID:
	    HYPRE_EuclidSolve(solver, parcsr_A, par_b, par_x);
	    break;
	  }
      }
    else
      {
	switch(type_preconditioner)
	  {
	  case BOOMER_AMG:
	    HYPRE_BoomerAMGSolveT(solver, parcsr_A, par_b, par_x);
	    break;
	  default :
	    cout << "Transpose preconditioning not available" << endl;
	    abort();
	  }
      }

    // then we expand the values in vector z (including shared rows)
    z.Zero();
    for (int i = 0; i < num.GetM(); i++)
      z(num(i)) = x_data[i];

    int nb_proc; MPI_Comm_size(comm, &nb_proc);    
    if (nb_proc > 1)
      AssembleVector(z, MPI_SUM, *ProcNumber, *DofNumber, comm, nodl_scalar, nb_u, 11);
  }

}

#define SELDON_FILE_HYPRE_CXX
#endif

