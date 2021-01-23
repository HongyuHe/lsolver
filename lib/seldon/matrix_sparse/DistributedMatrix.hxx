// Copyright (C) 2013-2015 INRIA
// Author(s): Marc Duruflé
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.

#ifndef SELDON_FILE_DISTRIBUTED_MATRIX_HXX

namespace Seldon
{

  //! class storing arrays needed for a distributed matrix
  class DistributedMatrixIntegerArray
  {
  public:
#ifdef SELDON_WITH_MPI
    //! row numbers shared with other processors
    /*!
      It is the same array required in a distributed vector.
      Some rows can be shared by several processors.
    */
    IVect OverlapRowNumbers;
    
    //! processor where each shared row should be assembled
    /*!
      For rows which have already been counted in other processors,
      you specify on which processor these rows are assembled
    */
    IVect OverlapProcNumbers;
    
    //! global row numbers 
    IVect GlobalRowNumbers;
    
    //! list of processors sharing rows with the current one
    IVect ProcSharingRows;
    
    //! for each processor sharing rows, list of "local numbers" shared
    /*!
      It is assumed that these local numbers are "sorted", such that
      local numbers shared with processor j on processor i are corresponding
      with local numbers shared with processor i on processor j
    */
    Vector<IVect> SharingRowNumbers;
    
    //! number of "scalar" unknowns
    int nodl_scalar, nb_unknowns_scal;
    
    //! total number of rows (on all processors)
    int nglob;

    //! local number of rows
    int nloc;
    
    //! MPI communicator
    MPI_Comm comm;

    void Clear();
    void SetData(int nl, int ng, int nodl, int nb_u, const MPI_Comm& comm_,
		 const IVect& overlap_row, const IVect& overlap_proc,
		 const IVect& global_rows, const IVect& proc_rows, const Vector<IVect>& sharing_rows);
    
    void Nullify();
#endif
    
  };


  //! Base class for distributed matrix over all the processors
  /*!
    In this class, distant non-zero entries are stored.
    It can contain non-zero entries corresponding to distant rows
    (i.e. the column belongs to the current processor, but
    the row belongs to another processor)
    or distant columns
    (i.e. the row belongs to the current processor, but
    the column belongs to another processor)

    Only one storage is implemented (equivalent of ArrayRowSparse).
   */
  template<class T>
  class DistributedMatrix_Base
  {
    template<class T0>
    friend class DistributedMatrix_Base;

  protected :
    //! row numbers shared with other processors
    /*!
      It is the same array required in a distributed vector.
      Some rows can be shared by several processors.
    */
    IVect* OverlapRowNumbers;
    
    //! processor where each shared row should be assembled
    /*!
      For rows which have already been counted in other processors,
      you specify on which processor these rows are assembled
    */
    IVect* OverlapProcNumbers;
    
    //! global row numbers 
    IVect* GlobalRowNumbers;
    
    //! list of processors sharing rows with the current one
    IVect* ProcSharingRows;
    
    //! for each processor sharing rows, list of "local numbers" shared
    /*!
      It is assumed that these local numbers are "sorted", such that
      local numbers shared with processor j on processor i are corresponding
      with local numbers shared with processor i on processor j
    */
    Vector<IVect>* SharingRowNumbers;
    
    //! number of "scalar" unknowns
    int nodl_scalar_, nb_unknowns_scal_;
    
    //! total number of rows (on all processors)
    int nglob_;
    
    //! MPI communicator
    MPI_Comm comm_;

    //! additional values on rows with non-local columns
    Vector<Vector<T, VectSparse>, VectFull,
           NewAlloc<Vector<T, VectSparse> > > dist_col;
    
    //! additional values on columns with non-local rows
    Vector<Vector<T, VectSparse>, VectFull,
                  NewAlloc<Vector<T, VectSparse> > > dist_row;
    
    //! distant processor for additional values
    Vector<IVect> proc_col, proc_row;
    
    //! global row/col numbers (needed for MltAdd)
    IVect global_row_to_recv, global_col_to_recv;
    IVect ptr_global_row_to_recv, ptr_global_col_to_recv;
    
    //! local row/col numbers (needed for MltAdd)
    Vector<IVect> local_row_to_send, local_col_to_send;
    
    //! processor numbers (needed for MltAdd)
    IVect proc_col_to_recv, proc_col_to_send,
      proc_row_to_recv, proc_row_to_send;
    
    //! if true local numbers are present in dist_row/dist_col 
    //! instead of global numbers
    bool local_number_distant_values;
    
    //! number of distant non-zero entries
    long size_max_distant_row, size_max_distant_col;
    
    // internal functions
    void EraseArrayForMltAdd();
    void SwitchToGlobalNumbers();
    
    template<class TypeDist>
    void SortAndAssembleDistantInteractions(TypeDist& dist_val,
					    Vector<IVect>& dist_proc,
                                            IVect& glob_num,
					    IVect& ptr_glob_num,
					    IVect& proc_glob,
                                            Vector<IVect>& local_num,
					    IVect& proc_local);
    
    template<class T2>
    void ScatterValues(const Vector<T2>& X, const IVect& num_recv,
                       const IVect&, const IVect& proc_recv,
                       const Vector<IVect>& num_send,
                       const IVect& proc_send, Vector<T2>& Xcol) const;
    
    template<class T2>
    void AssembleValues(const Vector<T2>& Xcol, const IVect& num_recv,
                        const IVect&, const IVect& proc_recv,
                        const Vector<IVect>& num_send,
                        const IVect& proc_send, Vector<T2>& X) const;
    
    void AssembleValuesMin(const IVect& Xcol, const IVect& Xcol_proc,
                           const IVect& num_recv, const IVect& ptr_num_recv,
                           const IVect& proc_recv,
                           const Vector<IVect>& num_send,
			   const IVect& proc_send,
                           IVect& Y, IVect& Yproc) const;
    
    void AssembleVecMin(Vector<int>& X, Vector<int>& Xproc) const;
    
    template<class T0, class TypeDist>
    void RemoveSmallEntryDistant(const T0&, TypeDist&, Vector<IVect>&);
    
    template<class T0> void GetRowSumDistantCol(Vector<T0>& vec_sum) const;
    template<class T0> void GetRowSumDistantRow(Vector<T0>& vec_sum) const;

    template<class T0> void GetColSumDistantCol(Vector<T0>& vec_sum) const;
    template<class T0> void GetColSumDistantRow(Vector<T0>& vec_sum) const;

    // internal functions for matrix-vector product    
    void PrepareMltAdd();
    
    template<class T2>
    void ScatterRowValues(const Vector<T2>& X, Vector<T2>& Xcol) const;

    template<class T2>
    void ScatterColValues(const Vector<T2>& X, Vector<T2>& Xcol) const;
    
    template<class T2>
    void AssembleRowValues(const Vector<T2>& Xrow, Vector<T2>& X) const;

    template<class T2>
    void AssembleColValues(const Vector<T2>& Xrow, Vector<T2>& X) const;

  public:
    template<class T2>
    void AssembleVec(Vector<T2>&) const;

  protected:
    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddCol(const SeldonTranspose& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;

    template<class T2, class Storage2, class Allocator2,
             class T4, class Storage4, class Allocator4>
    void MltAddRow(const SeldonTranspose& Trans,
                   const Vector<T2, Storage2, Allocator2>& X,
                   Vector<T4, Storage4, Allocator4>& Y) const;
    
    // internal static functions
    static void AddDistantValue(Vector<T, VectSparse>& dist_col_,
                                IVect& proc_col_,
                                int jglob, int proc2, const T& val);

    static void SendAndReceiveDistributed(const MPI_Comm& comm, IVect& nsend_int, Vector<IVect>& EntierToSend, 
                                          Vector<Vector<T> >& FloatToSend, IVect& nrecv_int,
                                          Vector<IVect>& EntierToRecv, Vector<Vector<T> >& FloatToRecv);
    
    
    static void AddReceivedInteractions(const MPI_Comm& comm, Matrix<T, General, ArrayRowSparse>& B,
                                        Vector<IVect>& EntierToRecv, Vector<Vector<T> >& FloatToRecv,
                                        IVect& nrecv_int, Vector<IVect>& EntierToSend,
                                        Vector<Vector<T> >& FloatToSend, IVect& nsend_int,
                                        IVect& Glob_to_local, const IVect& OverlappedCol,
                                        const IVect& OverlapProcNumber,
                                        Vector<IVect>& procB, bool reorder);
    
    static void AddReceivedInteractions(const MPI_Comm& comm, Matrix<T, General, ArrayColSparse>& B,
                                        Vector<IVect>& EntierToRecv, Vector<Vector<T> >& FloatToRecv,
                                        IVect& nrecv_int, Vector<IVect>& EntierToSend,
                                        Vector<Vector<T> >& FloatToSend, IVect& nsend_int,
                                        IVect& Glob_to_local, const IVect& OverlappedCol,
                                        const IVect& OverlapProcNumber,
                                        Vector<IVect>& procB, bool reorder);
    
    template<class TypeDist>
    static void EraseDistantEntries(MPI_Comm& comm, const Vector<bool>& IsRowDropped,
                                    const Vector<bool>& IsRowDroppedDistant,
                                    TypeDist& dist_row_, Vector<IVect>& proc_row_,
                                    TypeDist& dist_col_, Vector<IVect>& proc_col_);
    
  public :
    // constructors
    DistributedMatrix_Base();
    explicit DistributedMatrix_Base(int m, int n);
    
    // Inline methods
    MPI_Comm& GetCommunicator();
    const MPI_Comm& GetCommunicator() const;
  
    int GetLocalM() const;
    int GetLocalN() const;
    int GetGlobalM() const;
    int GetNodlScalar() const;
    int GetNbScalarUnknowns() const;

    void AddDistantInteraction(int i, int jglob, int proc,
                               const T& val);
    
    void AddRowDistantInteraction(int iglob, int j, int proc,
                                  const T& val);

    long GetMaxDataSizeDistantCol() const;
    long GetMaxDataSizeDistantRow() const;
    bool IsReadyForMltAdd() const;
    void SetReadyForMltAdd(bool all_zero = true);
    
    int GetDistantColSize(int i) const;
    int IndexGlobalCol(int i, int j) const;
    int ProcessorDistantCol(int i, int j) const;
    const T& ValueDistantCol(int i, int j) const;
    
    int GetDistantRowSize(int i) const;
    int IndexGlobalRow(int i, int j) const;
    int ProcessorDistantRow(int i, int j) const;
    const T& ValueDistantRow(int i, int j) const;
    
    // initialisation of pointers
    void Init(int n, IVect*, IVect*, IVect*,
              int, int, IVect*, Vector<IVect>*, const MPI_Comm&);

    void Init(DistributedMatrixIntegerArray&);
    
    template<class T0>
    void Init(const DistributedMatrix_Base<T0>&);
    
    void Init(Vector<IVect>&, IVect&, IVect&, IVect&,
              IVect&, Vector<IVect>&, const MPI_Comm&,
              bool distribute_row = true);

    void Init(IVect&, IVect&, IVect&,
              IVect&, Vector<IVect>&, const MPI_Comm&);
    
    // memory management
    void ReallocateDist(int m, int n);
    void Resize(int m, int n);
    void Clear();
    
    // basic functions
    DistributedMatrix_Base<T>&
    operator=(const DistributedMatrix_Base<T>& X);

    template<class T2>
    void Copy(const DistributedMatrix_Base<T2>& X);
    
    template<class T0>
    DistributedMatrix_Base<T>& operator *=(const T0& x);

    IVect& GetGlobalRowNumber();
    const IVect& GetGlobalRowNumber() const;
    IVect& GetOverlapRowNumber();
    const IVect& GetOverlapRowNumber() const;
    IVect& GetOverlapProcNumber();
    const IVect& GetOverlapProcNumber() const;
    IVect& GetProcessorSharingRows();
    const IVect& GetProcessorSharingRows() const;
    Vector<IVect>& GetSharingRowNumbers();
    const Vector<IVect>& GetSharingRowNumbers() const;

    size_t GetMemorySize() const;
    
    // convenient functions overloaded    
    long GetNonZeros() const;
    long GetDataSize() const;    
    
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);
    
    void SetIdentity();
    void Zero();
    void Fill();
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    
    void WriteText(ostream& FileStream, Vector<int>& Indow, Vector<int>& IndCol,
                   Vector<T>& Value, bool cplx) const;
    
    // functions called by matrix-vector products
    template<class T2, class T3, class T4, class Storage4, class Allocator4>
    void InitMltAdd(bool& proceed_distant_row, bool& proceed_distant_col,
                    const Vector<T2>& X, Vector<T2>& Xcol,
                    const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
                    Vector<T4, Storage4, Allocator4>& Yres) const;
    
    template<class T2, class T3, class T4, class Storage4, class Allocator4>
    void FinalizeMltAdd(bool proceed_distant_row, bool proceed_distant_col,
                        const Vector<T2>& X, Vector<T2>& Xcol, const T3& alpha,
                        const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
                        Vector<T4, Storage4, Allocator4>& Yres, bool assemble) const;
    
    template<class T2, class T3, class T4, class Storage4, class Allocator4>
    void InitMltAdd(bool& proceed_distant_row, bool& proceed_distant_col,
                    const SeldonTranspose& trans, const Vector<T2>& X, Vector<T2>& Xrow,
                    const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
                    Vector<T4, Storage4, Allocator4>& Yres) const;

    template<class T2, class T3, class T4, class Storage4, class Allocator4>
    void FinalizeMltAdd(bool proceed_distant_row, bool proceed_distant_col,
                        const SeldonTranspose& trans, const Vector<T2>& X, Vector<T2>& Xrow,
                        const T3& alpha, const T3& beta, Vector<T4, Storage4, Allocator4>& Y,
                        Vector<T4, Storage4, Allocator4>& Yres, bool assemble) const;
    
    void InitMltMin(Vector<int>& Y, Vector<int>& Yproc,
                    Vector<int>& Xcol, Vector<int>& Xcol_proc) const;

    void FinalizeMltMin(Vector<int>& Y, Vector<int>& Yproc,
                        Vector<int>& Xcol, Vector<int>& Xcol_proc) const;

    // methods called for various functions on matrices
    template<class T0, class T1>
    void AddDistributedMatrix(const T0& alpha,
                              const DistributedMatrix_Base<T1>& A);

    void GetMaxAbsDistant(typename ClassComplexType<T>::Treal& res) const;

    template<class T0>
    void AddRowSumDistant(Vector<T0>& vec_sum) const;

    template<class T0>
    void AddColSumDistant(Vector<T0>& vec_sum) const;

    template<class T0>
    void AddRowColSumDistant(Vector<T0>& sum_row, Vector<T0>& sum_col) const;
    
    void ExchangeParallelData(int& smax_row, int& smax_col, bool& local_number,
                              Vector<Vector<T, VectSparse>, VectFull,
                              NewAlloc<Vector<T, VectSparse> > >& dist_row_,
                              Vector<Vector<T, VectSparse>, VectFull,
                              NewAlloc<Vector<T, VectSparse> > >& dist_col_,
                              Vector<IVect>& proc_row_, Vector<IVect>& proc_col_,
                              IVect& global_row_to_recv_, IVect& global_col_to_recv_,
                              IVect& ptr_global_row_to_recv_, IVect& ptr_global_col_to_recv_,
                              Vector<IVect>& local_row_to_send_, Vector<IVect>& local_col_to_send_,
                              IVect& proc_row_to_recv_, IVect& proc_col_to_recv_,
                              IVect& proc_row_to_send_, IVect& proc_col_to_send_);

    void ConjugateDistant();
    void TransposeDistant(const DistributedMatrix_Base<T>& A);
    
    template<class T0>
    void ScaleLeftDistant(const Vector<T0>& Drow);
    
    template<class T0>
    void ScaleRightDistant(const Vector<T0>& Dcol);
    
    template<class T0, class T1>
    void ScaleDistant(const Vector<T0>& Drow, const Vector<T1>& Dcol);
    
    void EraseColDistant(const IVect& num, bool sym);
    void EraseRowDistant(const IVect& num, bool sym);

    template<class T1>
    void CopySubDistant(const DistributedMatrix_Base<T1>& A,
                        const IVect& row, const IVect& col, bool sym);
      

    // functions for assembling matrix
    void AssembleParallel(Matrix<T, General, ArrayRowSparse>& B, Vector<IVect>& procB,
                          Symmetric& sym, IVect& row_numbers, IVect& local_row_numbers,
                          IVect& OverlappedCol, bool sym_pattern, bool reorder);

    template<class Tint0, class Tint1>    
    static void ConvertToCSR(Matrix<T, General, ArrayRowSparse>& B, IVect& OverlappedCol,
                             Vector<Tint0>& PtrA, Vector<Tint1>& IndA, Vector<T>& ValA);

    void AssembleParallel(Matrix<T, General, ArrayColSparse>& B, Vector<IVect>& procB,
                          General& prop, IVect& col_numbers, IVect& local_col_numbers,
                          IVect& OverlappedCol, bool sym_pattern, bool reorder);

    template<class Tint0, class Tint1>    
    static void ConvertToCSC(Matrix<T, General, ArrayColSparse>& B, IVect& OverlappedCol,
                             Vector<Tint0>& PtrA, Vector<Tint1>& IndA, Vector<T>& ValA);
    
    template<class T0, class Allocator0>
    void GetDistributedRows(Matrix<T0, General,
			    ArrayRowSparse, Allocator0>& rows,
                            Vector<IVect>& proc, bool sym) const;
    
    template<class T0, class Allocator0>
    void GetDistributedColumns(Matrix<T0, General, ArrayColSparse, Allocator0>& B,
                               Vector<IVect>& procB, Vector<long>& Ptr, IVect& Ind,
                               Vector<T0>& Val, bool sym_pattern) const;
    
  };


  //! matrix distributed over all the processors
  /*!
    The local part of the matrix (with local rows and columns)
    is stored in the parent class Matrix<T, Prop, Storage, Allocator>
    This part can be stored at the user convenience (RowSparse, ArrayRowSymSparse, etc)
    The distant part of the matrix (with distant rows or distant columns)
    is stored in the parent class DistributedMatrix_Base

    Most of the functions that apply to standard matrices
    are overloaded for the class DistributedMatrix such that
    it should work correctly in parallel.
   */
  template<class T, class Prop, class Storage, class Allocator
	   = typename SeldonDefaultAllocator<Storage, T>::allocator>
  class DistributedMatrix : public Matrix<T, Prop, Storage, Allocator>,
                            public DistributedMatrix_Base<T>
  {
    template<class T0, class Prop0, class Storage0, class Allocator0>
    friend class DistributedMatrix;
    
  public :
    // constructors
    DistributedMatrix();
    explicit DistributedMatrix(int m, int n);

    // Inline methods
    void AddDistantInteraction(int i, int jglob, int proc,
                               const T& val);
    
    void AddRowDistantInteraction(int iglob, int j, int proc,
                                  const T& val);

    // memory management
    void Reallocate(int m, int n);
    void Resize(int m, int n);
    void Clear();
    void ClearLocal();
    
    // basic functions
    template<class T2, class Prop2, class Storage2, class Allocator2>
    DistributedMatrix<T, Prop, Storage, Allocator>&
    operator=(const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& X);
    
    template<class T0>
    DistributedMatrix<T, Prop, Storage, Allocator>& operator *=(const T0& x);

    size_t GetMemorySize() const;
    
    // convenient functions
    long GetNonZeros() const;
    long GetDataSize() const;    
    
    template<class T0>
    void RemoveSmallEntry(const T0& epsilon);
    
    void SetIdentity();
    void Zero();
    void Fill();
    
    template<class T0>
    void Fill(const T0& x);
    
    void FillRand();
    
    void Write(string FileName) const;
    void Write(ostream& FileStream) const;
    void WriteText(string FileName, bool cplx = false) const;
    void WriteText(ostream& FileStream, bool cplx = false) const;
    void Read(string FileName);
    void Read(istream& FileStream);
    void ReadText(string FileName, bool cplx = false);
    void ReadText(istream& FileStream, bool cplx = false);
    
    // functions for assembling matrix
    template<class T0, class Allocator0>
    void GetDistributedRows(Matrix<T0, General,
			    ArrayRowSparse, Allocator0>& rows,
                            Vector<IVect>& proc) const;
    
    template<class T0, class Allocator0>
    void GetDistributedColumns(Matrix<T0, General,
			       ArrayColSparse, Allocator0>& rows,
			       Vector<IVect>&, bool sym_pattern) const;

#ifdef SELDON_WITH_VIRTUAL
    // virtual inline methods
    typedef typename ClassComplexType<T>::Treal Treal;
    typedef typename ClassComplexType<T>::Tcplx Tcplx;
    
    virtual void ApplySor(const SeldonTranspose&, Vector<Treal>& x, const Vector<Treal>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;
    
    virtual void ApplySor(const SeldonTranspose&, Vector<Tcplx>& x, const Vector<Tcplx>& r,
			  const typename ClassComplexType<T>::Treal& omega,
			  int nb_iter, int stage_ssor) const;
    
    virtual void MltAddVector(const Treal& alpha, const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;

    virtual void MltAddVector(const Treal& alpha, const SeldonTranspose&,
			      const Vector<Treal>& x,
			      const Treal& beta, Vector<Treal>& y) const;

    virtual void MltAddVector(const Tcplx& alpha, const SeldonTranspose&,
			      const Vector<Tcplx>& x,
			      const Tcplx& beta, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const Vector<Treal>& x, Vector<Treal>& y) const;
    virtual void MltVector(const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
    
    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Treal>& x, Vector<Treal>& y) const;

    virtual void MltVector(const SeldonTranspose&,
			   const Vector<Tcplx>& x, Vector<Tcplx>& y) const;
#endif
    
  };


  // matrix vector products  
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Storage4, Allocator4>& Yres, bool assemble = true);
  
  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Storage2, class Allocator2, class T3,
           class T4, class Storage4, class Allocator4>
  void MltAddVector(const T0& alpha,
		    const SeldonTranspose& Trans,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		    const Vector<T2, Storage2, Allocator2>& X,
		    const T3& beta,
		    Vector<T4, Storage4, Allocator4>& Y, bool assemble = true);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void MltVector(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& M,
		 const Vector<T1, Storage1, Allocator1>& X,
		 Vector<T2, Storage2, Allocator2>& Y, bool assemble = true);
  
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void MltVector(const T3& alpha,
		 const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);
  
  template <class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void MltVector(const SeldonTranspose& Trans,
		 const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
		 const Vector<T2, Storage2, Allocator2>& X,
		 Vector<T3, Storage3, Allocator3>& Y, bool assemble = true);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1>
  void MltScalar(const T0& alpha,
		 DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc);

  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc);

#ifdef SELDON_FILE_MATRIX_ARRAY_COMPLEX_SPARSE_HXX
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc);
  
  template<class T1, class Prop1, class Allocator1>
  void MltMin(const Matrix<T1, Prop1,
	      ArrayRowSymComplexSparse, Allocator1>& A,
              const IVect& global, IVect& Y, IVect& Yproc);
#endif
  
  template<class T1, class Prop1, class Storage1, class Allocator1>
  void MltMin(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& M,
              IVect& Y, IVect& Yproc);

  // matrix functions

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void AddMatrix(const T0& alpha,
		 const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A, 
		 DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  MaxAbs(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetRowSum(Vector<T0>& vec_sum,
                 const DistributedMatrix<T, Prop, Storage, Allocator>& A);

  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetColSum(Vector<T0>& vec_sum,
                 const DistributedMatrix<T, Prop, Storage, Allocator>& A);

  template<class T0, class T, class Prop, class Storage, class Allocator>
  void GetRowColSum(Vector<T0>& sum_row,
                    Vector<T0>& sum_col,
                    const DistributedMatrix<T, Prop, Storage, Allocator> & A);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  Norm1(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  typename ClassComplexType<T1>::Treal
  NormInf(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A); 

  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Transpose(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Conjugate(DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  void Transpose(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
                 DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B);

  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(const DistributedMatrix<T, Prop,
		     Storage, Allocator>& A,
                     DistributedMatrix<T, Prop, Storage, Allocator>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(DistributedMatrix<T, Prop, Storage, Allocator>& A);

  // matrix-matrix products
  
  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltMatrix(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		 const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
		 DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C);

  template<class T0,
           class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2,
           class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAddMatrix(const T0& alpha,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		    const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
		    const T3& beta,
		    DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C);

  template<class T0, class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2, class T3,
           class T4, class Prop4, class Storage4, class Allocator4>
  void MltAddMatrix(const T0& alpha, const SeldonTranspose& transA,
		    const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
		    const SeldonTranspose& transB,
		    const DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B,
		    const T3& beta,
		    DistributedMatrix<T4, Prop4, Storage4, Allocator4>& C);

  // Matrix functions
  
  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetRow(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void GetCol(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
              int i, Vector<T1, VectSparse, Allocator1>& X);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetRow(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Allocator1>
  void SetCol(const Vector<T1, VectSparse, Allocator1>& X,
              int i, DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void ApplyPermutation(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                        const Vector<int>& row_perm,
                        const Vector<int>& col_perm);

  template<class T, class Prop, class Storage, class Allocator>
  void 
  ApplyInversePermutation(DistributedMatrix<T, Prop, Storage, Allocator>& A,
			  const Vector<int>& row_perm,
			  const Vector<int>& col_perm);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 2);

  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2, class T3>
  void SorVector(const SeldonTranspose& transM,
		 const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		 Vector<T2, Storage2, Allocator2>& X,
		 const Vector<T1, Storage1, Allocator1>& B,
		 const T3& omega, int iter, int type_ssor = 3);

  template<class T1, class Prop, class Storage, class Allocator,
	   class T2, class Allocator2, class Allocator3>
  void GetCol(const DistributedMatrix<T1, Prop, Storage, Allocator>& A,
              const IVect& col_number,
	      Vector<Vector<T2, VectSparse, Allocator2>,
	      VectSparse, Allocator3>& V);

  template<class T, class Prop1, class Storage1, class Allocator1,
           class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<T, Prop2, Storage2, Allocator2>& B);

  template<class T, class Prop1, class Storage1, class Allocator1,
           class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<T, Prop2, Storage2, Allocator2>& B);
  
  template<class T, class Prop1, class Storage1, class Allocator1,
           class Prop2, class Storage2, class Allocator2>
  void Copy(const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
            DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B);

  template<class T1, class Prop1, class Storage1, class Allocator1,
           class T2, class Prop2, class Storage2, class Allocator2>
  void
  CopyReal(const DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A,
	   DistributedMatrix<T2, Prop2, Storage2, Allocator2>& B);
  
  template<class T, class Prop, class Storage, class Allocator>
  typename ClassComplexType<T>::Treal
  NormFro(const DistributedMatrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Storage, class Allocator,
           class T1, class Allocator1>
  void ScaleLeftMatrix(DistributedMatrix<T, General, Storage, Allocator>& A,
                       const Vector<T1, VectFull, Allocator1>& Drow);

  template<class T, class Storage, class Allocator,
           class T1, class Allocator1>
  void ScaleRightMatrix(DistributedMatrix<T, General, Storage, Allocator>& A,
                        const Vector<T1, VectFull, Allocator1>& Dcol);

  template<class T, class Prop, class Storage, class Allocator,
           class T1, class Allocator1, class T2, class Allocator2>
  void ScaleMatrix(DistributedMatrix<T, Prop, Storage, Allocator>& A,
                   const Vector<T1, VectFull, Allocator1>& Drow,
                   const Vector<T2, VectFull, Allocator2>& Dcol);

  template<class Prop, class Storage, class Alloc, class Tint0, class Tint1, class T>
  void AssembleDistributed(DistributedMatrix<T, Prop, Storage, Alloc>& A,
			   Symmetric& sym, const MPI_Comm& comm,
                           IVect& row_numbers, IVect& local_row_numbers,
			   Vector<Tint0>& PtrA, Vector<Tint1>& IndA,
                           Vector<T>& ValA, bool sym_pattern, bool reorder = false);

  template<class Prop, class Storage, class Alloc, class Tint0, class Tint1, class T>
  void AssembleDistributed(DistributedMatrix<T, Prop, Storage, Alloc>& A,
			   General& prop, const MPI_Comm& comm,
                           IVect& col_numbers, IVect& local_col_numbers,
                           Vector<Tint0>& PtrA, Vector<Tint1>& IndA,
                           Vector<T>& ValA, bool sym_pattern, bool reorder = false);

  template<class TypeDist>
  void EraseDistantEntries(MPI_Comm& comm, const Vector<bool>& IsRowDropped,
                           const Vector<bool>& IsRowDroppedDistant,
                           TypeDist& dist_row, Vector<IVect>& proc_row,
                           TypeDist& dist_col, Vector<IVect>& proc_col);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  void EraseCol(const IVect& num,
                DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T1, class Prop1, class Storage1, class Allocator1>
  void EraseRow(const IVect& num,
                DistributedMatrix<T1, Prop1, Storage1, Allocator1>& A);

  template<class T0, class Prop0, class Storage0, class Allocator0,
           class T1, class Prop1, class Storage1, class Allocator1>
  void
  CopySubMatrix(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
		const IVect& row, const IVect& col,
                DistributedMatrix<T1, Prop1, Storage1, Allocator1>& B);

  /****************************************
   * Mlt, MltAdd for distributed matrices *
   ****************************************/

  // functions present in DistributedMatrixInline.cxx
  template<class T, class Prop, class Storage, class Allocator>
  void Mlt(const T& alpha,
           DistributedMatrix<T, Prop, Storage, Allocator>& A);

  template<class T, class Prop, class Storage, class Allocator>
  void Mlt(const T& alpha,
           DistributedMatrix<complex<T>, Prop, Storage, Allocator>& A);

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const T& alpha,
           const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
           DistributedMatrix<T, Prop2, Storage2, Allocator2>& B);
  
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const complex<T>& alpha,
           const DistributedMatrix<T, Prop1, Storage1, Allocator1>& A,
           DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B);

  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const T& alpha,
           const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
           DistributedMatrix<complex<T>, Prop2, Storage2, Allocator2>& B);
  
  template <class T,
	    class Prop1, class Storage1, class Allocator1,
	    class Prop2, class Storage2, class Allocator2>
  void Add(const T& alpha,
           const DistributedMatrix<complex<T>, Prop1, Storage1, Allocator1>& A,
           DistributedMatrix<T, Prop2, Storage2, Allocator2>& B);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>&,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<T0, Prop0, Storage0, Allocator0>&,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y,
	   bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>&,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Vector<complex<T0> >& X, Vector<complex<T0> >& Y,
	   bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void Mlt(const SeldonTranspose& trans,
	   const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	   const Vector<T0>& X, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y,
	      bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const complex<T0>& alpha,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const complex<T0>& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<T0, Prop0, Storage0, Allocator0>& A,
	      const Vector<complex<T0> >& X, const complex<T0>& beta,
	      Vector<complex<T0> >& Y, bool assemble = true);
  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void MltAdd(const T0& alpha, const SeldonTranspose& trans,
	      const DistributedMatrix<complex<T0>, Prop0, Storage0, Allocator0>& A,
	      const Vector<T0>& X, const T0& beta, Vector<T0>& Y, bool assemble = true);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const T& omega, int iter, int type_ssor = 2);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const SeldonTranspose&,
	   const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const T& omega, int iter, int type_ssor = 3);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<complex<T>, Storage2, Allocator2>& Y,
	   const Vector<complex<T>, Storage1, Allocator1>& X,
	   const T& omega, int iter, int type_ssor = 2);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const SeldonTranspose&,
	   const DistributedMatrix<T, Prop0, Storage0, Allocator0>& M,
	   Vector<complex<T>, Storage2, Allocator2>& Y,
	   const Vector<complex<T>, Storage1, Allocator1>& X,
	   const T& omega, int iter, int type_ssor = 3);
  
  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const DistributedMatrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const T& omega, int iter, int type_ssor = 2);

  template <class T, class Prop0, class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>	    
  void SOR(const SeldonTranspose&,
	   const DistributedMatrix<complex<T>, Prop0, Storage0, Allocator0>& M,
	   Vector<T, Storage2, Allocator2>& Y,
	   const Vector<T, Storage1, Allocator1>& X,
	   const T& omega, int iter, int type_ssor = 3);

}

#define SELDON_FILE_DISTRIBUTED_MATRIX_HXX
#endif

