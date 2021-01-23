#ifndef SELDON_FILE_HYPRE_HXX

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"

namespace Seldon
{

  //! class interfacing Hypre with Seldon (as a preconditioner)
  /*!
    This class is declared template but can be used only for T = HYPRE_Complex
   */
  template<class T>
  class HyprePreconditioner : public Preconditioner_Base<T>
  {
  protected:
    //! distributed matrix stored as in hypre
    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    
    //! hypre solver 
    HYPRE_Solver solver;
    bool solver_created;
    
    //! intermediary vectors x and b
    HYPRE_IJVector vec_b, vec_x;
    HYPRE_ParVector par_x;
    HYPRE_ParVector par_b;
    
    //! row numbers of original rows (that belong to only one processor)
    //! In Seldon row numbers can be shared, so we have to make the conversion
    //! between seldon vectors to "hypre" vectors (where the rows are not shared) through this array
    Vector<int> local_rows;
    
    //! list of processors implied in the assembling step
    Vector<int>* ProcNumber;

    //! list of rows implied in the assembling step
    Vector<Vector<int> >* DofNumber;
    
    //! number of rows for an unknown and the number of unknowns
    int nodl_scalar, nb_u;
    
    //! MPI communicator
    MPI_Comm comm;
    
    //! verbose level
    int print_level;
    
    //! hypre preconditioning selected
    int type_preconditioner;
    
    // parameters for Boomer AMG
    int amg_max_levels, amg_num_sweeps, amg_smoother;  
    
    // parameters for ParaSails
    HYPRE_Real sai_filter, sai_threshold;
    int sai_max_levels, sai_sym;
    
    // parameters for Euclid
    int euclid_level;
    bool euclid_use_ilut;
    HYPRE_Real euclid_threshold, euclid_droptol;
    
  public:
    
    //! available interfaced preconditionings in Hypre (for CSR matrices)
    enum {BOOMER_AMG, PARASAILS, EUCLID, AMS};
    
    //! available smoothers
    enum {JACOBI=0, GS_SEQ=1, GS_PAR_SEQ=2, HYBRID_GS_BACKWARD=3, HYBRID_GS_FORWARD=4,
	  HYBRID_GS_SYMMETRIC=6, L1_GAUSS_SEIDEL=8, CHEBYSHEV=16, FCF_JACOBI=17, L1_JACOBI=18};
    
    
    HyprePreconditioner();
    ~HyprePreconditioner();
    
    void Clear();

    // inline methods
    void SetPreconditioner(int type);
    void SetSmoother(int type);
    void SetLevelEuclid(int lvl);
    void ShowMessages();

    // specifying preconditioning with a line of a data file
    void SetInputPreconditioning(const string&, const Vector<string>&);
    
    // main method constructing the preconditioning
    template<class Prop, class Storage, class Allocator>
    void ConstructPreconditioner(DistributedMatrix<HYPRE_Complex, Prop, Storage, Allocator>& A0,
				 bool keep_matrix = false);

  protected:
    // internal method
    void FinalizePreconditioner(Vector<int>& row_numbers, Vector<int>& size_rows,
				Vector<int>& cols, Vector<HYPRE_Complex>& values);

  public :
    // methods to apply preconditioning
    void Solve(const SeldonTranspose& trans, const VirtualMatrix<T>& A,
	       const Vector<T>& r, Vector<T>& z);

    void Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z);
    void TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z);
    
  };
  
  
}

#define SELDON_FILE_HYPRE_HXX
#endif

