# main part of Seldon
OBJ_SELDON := $(CHESELDON)/lib/Compil/Seldon/Array3D$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixPointers$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixPacked$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixHermSymTriang$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixSparse$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Vector$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Vector2$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixVector$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrix$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixDense$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixArray$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixConversion$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/PermutationScalingMatrix$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/RelaxationMatrixVector$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/IOMatrixMarket$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/IlutPreconditioning$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/IterativeSolver$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/check_dim$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/TinyVector$(SELD_SUFFIX).o

# complex matrices
OBJ_SELDON := $(OBJ_SELDON) $(CHESELDON)/lib/Compil/Seldon/FunctionsMatVectComplex$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixComplex$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixComplexConversions$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/MatrixComplexSparse$(SELD_SUFFIX).o

# Blas functions
ifeq ($(USE_BLAS),YES)
  OBJ_BLAS := $(CHESELDON)/lib/Compil/Seldon/Blas$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o
else
  OBJ_BLAS := $(CHESELDON)/lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o
endif
OBJ_SELDON := $(OBJ_SELDON) $(OBJ_BLAS)
OBJ_BLAS := $(CHESELDON)/lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Blas$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o

# MPI functions
ifeq ($(USE_MPI),YES)
  OBJ_SELDON := $(OBJ_SELDON) $(CHESELDON)/lib/Compil/Seldon/DistributedVector$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/DistributedMatrix$(SELD_SUFFIX).o
endif

# interface with eigenvalue solvers
# Arpack is always needed in that case
OBJ_EIG := $(CHESELDON)/lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o
ifeq ($(USE_ARPACK),YES)
  OBJ_EIG := $(OBJ_EIG) $(CHESELDON)/lib/Compil/Seldon/Arpack$(SELD_SUFFIX).o
  ifeq ($(USE_FEAST),YES)
    OBJ_EIG := $(OBJ_EIG) $(CHESELDON)/lib/Compil/Seldon/Feast$(SELD_SUFFIX).o
  endif
  ifeq ($(USE_ANASAZI),YES)
    OBJ_EIG := $(OBJ_EIG) $(CHESELDON)/lib/Compil/Seldon/Anasazi$(SELD_SUFFIX).o
  endif
  ifeq ($(USE_SLEPC),YES)
    OBJ_EIG := $(OBJ_EIG) $(CHESELDON)/lib/Compil/Seldon/Slepc$(SELD_SUFFIX).o
  endif
endif

OBJ_SELDON := $(OBJ_SELDON) $(OBJ_EIG)

# interfaces with direct solvers
OBJ_SELDON_SOLVE := $(CHESELDON)/lib/Compil/Seldon/CholeskySolver$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/SparseSeldonSolver$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/SparseDirectSolver$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/DistributedSolver$(SELD_SUFFIX).o

ifeq ($(USE_MUMPS),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o
endif

ifeq ($(USE_PARDISO),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/Pardiso$(SELD_SUFFIX).o
endif

ifeq ($(USE_PASTIX),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o
endif

ifeq ($(USE_PASTIX),INT64)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o
endif

ifeq ($(USE_UMFPACK),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/UmfPack$(SELD_SUFFIX).o $(CHESELDON)/lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o
endif

ifeq ($(USE_SUPERLU),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/SuperLU$(SELD_SUFFIX).o
endif

ifeq ($(USE_CHOLMOD),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o
endif

ifeq ($(USE_WSMP),YES)
  OBJ_SELDON_SOLVE := $(OBJ_SELDON_SOLVE) $(CHESELDON)/lib/Compil/Seldon/Wsmp$(SELD_SUFFIX).o
endif

OBJ_SELDON := $(OBJ_SELDON) $(OBJ_SELDON_SOLVE)

# list of dependances for each file .cpp contained in $(CHESELDON)/lib/Compil/Seldon
$(CHESELDON)/lib/Compil/Seldon/Array3D$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Array3D.cpp $(CHESELDON)/array/Array3D.cxx $(CHESELDON)/array/Array4D.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/BandMatrix$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/BandMatrix.cpp $(CHESELDON)/matrix_sparse/BandMatrix.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Common.cpp $(CHESELDON)/share/Errors.cxx $(CHESELDON)/share/MatrixFlag.cxx $(CHESELDON)/share/Common.cxx $(CHESELDON)/share/Allocator.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/DistributedVector$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/DistributedVector.cpp $(CHESELDON)/share/MpiCommunication.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/DistributedMatrix$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/DistributedMatrix.cpp $(CHESELDON)/matrix_sparse/DistributedMatrix.cxx $(CHESELDON)/matrix_sparse/DistributedMatrixFunction.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixPointers$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixPointers.cpp $(CHESELDON)/matrix/Matrix_Pointers.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixPacked$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixPacked.cpp $(CHESELDON)/matrix/Matrix_HermPacked.cxx $(CHESELDON)/matrix/Matrix_SymPacked.cxx $(CHESELDON)/matrix/Matrix_TriangPacked.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixHermSymTriang$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixHermSymTriang.cpp $(CHESELDON)/matrix/Matrix_Hermitian.cxx $(CHESELDON)/matrix/Matrix_Symmetric.cxx $(CHESELDON)/matrix/Matrix_Triangular.cxx $(CHESELDON)/matrix/Matrix_Base.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixSparse$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixSparse.cpp $(CHESELDON)/matrix_sparse/Matrix_Sparse.cxx $(CHESELDON)/matrix_sparse/Matrix_SymSparse.cxx $(CHESELDON)/matrix_sparse/Matrix_ArraySparse.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Vector$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Vector.cpp $(CHESELDON)/vector/Vector.cxx $(CHESELDON)/vector/SparseVector.cxx $(CHESELDON)/vector/Functions_Arrays.cxx $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Vector2$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Vector2.cpp $(CHESELDON)/vector/Vector2.cxx $(CHESELDON)/vector/Vector3.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixVector$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixVector.cpp $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx $(CHESELDON)/computation/interfaces/Mkl_Sparse.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/FunctionsMatrix$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrix.cpp $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx $(CHESELDON)/computation/interfaces/Mkl_Sparse.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixDense$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixDense.cpp $(CHESELDON)/matrix/Functions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixArray$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixArray.cpp $(CHESELDON)/matrix_sparse/Functions_MatrixArray.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/FunctionsMatVectComplex$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/FunctionsMatVectComplex.cpp $(CHESELDON)/matrix_sparse/complex/Functions_MatVectComplex.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixComplex$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/FunctionsMatrixComplex.cpp $(CHESELDON)/matrix_sparse/complex/Functions_MatrixComplex.cxx $(CHESELDON)/matrix_sparse/Functions_MatrixArray.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixComplexConversions$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixComplexConversions.cpp $(CHESELDON)/matrix_sparse/complex/Matrix_ComplexConversions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixComplexSparse$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixComplexSparse.cpp $(CHESELDON)/matrix_sparse/complex/Matrix_ArrayComplexSparse.cxx $(CHESELDON)/matrix_sparse/complex/Matrix_ComplexSparse.cxx $(CHESELDON)/matrix_sparse/complex/Matrix_SymComplexSparse.cxx $(CHESELDON)/matrix_sparse/Matrix_Sparse.cxx $(CHESELDON)/matrix_sparse/IOMatrixMarket.cxx $(CHESELDON)/matrix/Matrix_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/MatrixConversion$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/MatrixConversion.cpp $(CHESELDON)/matrix_sparse/Matrix_Conversions.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/PermutationScalingMatrix$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/PermutationScalingMatrix.cpp $(CHESELDON)/matrix_sparse/Permutation_ScalingMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/RelaxationMatrixVector$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/RelaxationMatrixVector.cpp $(CHESELDON)/matrix_sparse/Relaxation_MatVect.cxx $(CHESELDON)/matrix_sparse/complex/Functions_MatVectComplex.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Blas$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Blas.cpp $(CHESELDON)/computation/interfaces/Blas_1.cxx $(CHESELDON)/computation/interfaces/Blas_2.cxx $(CHESELDON)/computation/interfaces/Blas_3.cxx $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/BlasMpfr$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/BlasMpfr.cpp $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx 
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/CholeskySolver$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/CholeskySolver.cpp $(CHESELDON)/computation/solver/SparseCholeskyFactorisation.cxx $(CHESELDON)/computation/solver/DistributedCholeskySolver.cxx 
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/EigenvalueSolver$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/EigenvalueSolver.cpp $(CHESELDON)/computation/interfaces/eigenvalue/VirtualEigenvalueSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Arpack$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Arpack.cpp $(CHESELDON)/computation/interfaces/eigenvalue/ArpackSolver.cxx $(CHESELDON)/computation/interfaces/eigenvalue/Arpack.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Feast$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Feast.cpp $(CHESELDON)/computation/interfaces/eigenvalue/Feast.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Anasazi$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Anasazi.cpp $(CHESELDON)/computation/interfaces/eigenvalue/Anasazi.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Slepc$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Slepc.cpp $(CHESELDON)/computation/interfaces/eigenvalue/Slepc.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/IlutPreconditioning$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/IlutPreconditioning.cpp $(CHESELDON)/computation/solver/preconditioner/IlutPreconditioning.cxx $(CHESELDON)/computation/solver/preconditioner/SymmetricIlutPreconditioning.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/IOMatrixMarket$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/IOMatrixMarket.cpp $(CHESELDON)/matrix_sparse/IOMatrixMarket.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/IterativeSolver$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/IterativeSolver.cpp $(CHESELDON)/computation/solver/iterative/Iterative.cxx $(CHESELDON)/computation/solver/preconditioner/Precond_Ssor.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Lapack$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Lapack.cpp $(CHESELDON)/computation/interfaces/Lapack_LinearEquations.cxx $(CHESELDON)/computation/interfaces/Lapack_LeastSquares.cxx $(CHESELDON)/computation/interfaces/Lapack_Eigenvalues.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/LapackMpfr$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/LapackMpfr.cpp $(CHESELDON)/computation/basic_functions/Functions_Base.cxx src/Algebra/FactorisationLU.cxx src/Algebra/Eigenvalue.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Cholmod$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Cholmod.cpp $(CHESELDON)/computation/interfaces/direct/Cholmod.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Mumps$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Mumps.cpp $(CHESELDON)/computation/interfaces/direct/Mumps.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Pastix$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Pastix.cpp $(CHESELDON)/computation/interfaces/direct/Pastix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Pardiso$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Pardiso.cpp $(CHESELDON)/computation/interfaces/direct/Pardiso.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/UmfPack$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/UmfPack.cpp $(CHESELDON)/computation/interfaces/direct/UmfPack.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/SuperLU$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/SuperLU.cpp $(CHESELDON)/computation/interfaces/direct/SuperLU.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/Wsmp$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/Wsmp.cpp $(CHESELDON)/computation/interfaces/direct/Wsmp.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/SparseSeldonSolver$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/SparseSeldonSolver.cpp $(CHESELDON)/computation/solver/Ordering.cxx $(CHESELDON)/computation/solver/SparseSolver.cxx $(CHESELDON)/computation/basic_functions/Functions_Base.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/SparseDirectSolver$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/SparseDirectSolver.cpp $(CHESELDON)/computation/interfaces/direct/SparseDirectSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/DistributedSolver$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/DistributedSolver.cpp $(CHESELDON)/computation/solver/DistributedSolver.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/check_dim$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/check_dim.cpp $(CHESELDON)/computation/basic_functions/Functions_Vector.cxx $(CHESELDON)/computation/basic_functions/Functions_MatVect.cxx $(CHESELDON)/computation/basic_functions/Functions_Matrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)

$(CHESELDON)/lib/Compil/Seldon/TinyVector$(SELD_SUFFIX).o : $(CHESELDON)/lib/Compil/Seldon/TinyVector.cpp $(CHESELDON)/vector/TinyVector.cxx $(CHESELDON)/matrix/TinyMatrix.cxx
	$(VERBOSE)$(CC) -c $(INCLUDE) $(FLAGS_OBJ) $< -o $@ $(AGRESSIVE_OPTIM)


$(LIB_SELDON_STATIC) :  $(OBJ_SELDON)
	ar rv $(LIB_SELDON_STATIC) $(OBJ_SELDON)

$(LIB_SELDON) :  $(OBJ_SELDON)
	$(CC) -shared -Wl,-soname,libseldon$(SELD_SUFFIX).so -rdynamic -o $(LIB_SELDON) $(OBJ_SELDON)

# target to remove all the object files
cleanlib :
	rm -f $(LIB_SELDON) $(LIB_SELDON_STATIC) $(OBJ_SELDON)

# target to remove files related to the interface with the direct solvers
cleansolve :
	rm -f $(OBJ_SELDON_SOLVE)

# target to remove files related to the interface with the eigenvalue solvers
cleaneig :
	rm -f $(OBJ_EIG) $(CHESELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o

# target to switch to SELDON_WITH_ABORT
cleanabort :
	rm -f $(CHESELDON)/lib/Compil/Seldon/Common$(SELD_SUFFIX).o

# Uncomment the following line in order to detect the file .o which gives bad result
#LIB_SELDON := $(OBJ_SELDON)

ifeq ($(STATIC_COMPILATION),YES)
  LIB_SELDON := $(LIB_SELDON_STATIC)
endif

ifeq ($(SEPARED_COMPIL),YES)
  LIB := $(LIB_SELDON) $(LIB)
endif
