#ifndef SELDON_FILE_HYPRE_INLINE_CXX

namespace Seldon
{

  //! sets which preconditioning to use (among choices proposed in hypre)
  template<class T>
  inline void HyprePreconditioner<T>::SetPreconditioner(int type)
  {
    type_preconditioner = type;
  }


  //! sets smoother to use for AMG
  template<class T>
  inline void HyprePreconditioner<T>::SetSmoother(int type)
  {
    amg_smoother = type;
  }


  //! sets the level k of ILU(k) factorization
  template<class T>
  inline void HyprePreconditioner<T>::SetLevelEuclid(int lvl)
  {
    euclid_level = lvl;
  }
  

  //! shows informations about preconditioning
  template<class T>
  inline void HyprePreconditioner<T>::ShowMessages()
  {
    print_level = 3;
  }
  

  //! applies preconditioning
  template<class T>
  inline void HyprePreconditioner<T>
  ::Solve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
  {
    Solve(SeldonNoTrans, A, r, z);
  }


  //! applies transpose preconditioning
  template<class T>
  inline void HyprePreconditioner<T>
  ::TransSolve(const VirtualMatrix<T>& A, const Vector<T>& r, Vector<T>& z)
  {
    Solve(SeldonTrans, A, r, z);
  }
  
}

#define SELDON_FILE_HYPRE_INLINE_CXX
#endif

