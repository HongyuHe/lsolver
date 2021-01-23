// Copyright (C) 2001-2009 Vivien Mallet
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


#include "SeldonLib.hxx"

using namespace Seldon;

typedef complex<double> cpx;
double epsilon(1e-12);


template<class T, class Prop, class Storage,
	 class Allocator1, class Allocator2, class Allocator3>
void Solve(Matrix<T, Prop, Storage, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x,
	   const Vector<T, VectFull, Allocator3>& b)
{
  SparseDirectSolver<T> mat_lu;
  
  // The initial matrix is erased during the factorization process
  // if you want to keep it, type mat_lu.Factorize(A, true);
  mat_lu.Factorize(A);
  x = b;
  mat_lu.Solve(x);
}

template<class T, class Prop, class Storage,
	 class Allocator1, class Allocator2, class Allocator3>
void Solve(Matrix<T, Prop, Storage, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x,
	   const Vector<T, VectFull, Allocator3>& b,
	   Vector<T, VectFull, Allocator3>& bt)
{
  SparseDirectSolver<T> mat_lu;
  
  // The initial matrix is erased during the factorization process
  mat_lu.Factorize(A);
  x = b;
  mat_lu.Solve(x);
  mat_lu.Solve(SeldonTrans, bt);
}

template<class Vector1>
bool CheckSolution(const Vector1& x)
{
  typename Vector1::value_type one;
  SetComplexOne(one);
  Vector1 xt(x.GetM());
  xt.Fill();
  Add(-one, x, xt);
  return Norm2(xt) <= epsilon;
}

template<class Vector1>
bool CheckSolution(const Vector1& x, const Vector1& x_ref, double& err)
{
  typename Vector1::value_type one;
  SetComplexOne(one);
  Vector1 xt(x);
  Add(-one, x_ref, xt);
  err = Norm2(xt)/Norm2(x_ref);
  return ( err <= epsilon );
}


int main(int argc, char **argv)
{
#ifdef SELDON_WITH_MPI
  MPI_Init(&argc, &argv);
#endif

  // Number of rows.
  int n = 5;

  // Results.
  bool success_rowsym_real, success_rowsym_complex, success_row_real, success_row_complex;
  
  // for transpose
  bool successT_row_real, successT_row_complex;

  // We test RowSymSparse (real numbers).
  {
    // construction of a matrix
    Matrix<double, Symmetric, RowSymSparse> A;
    int colind_[] = {0, 1, 1, 2, 2, 3, 3, 4, 4};
    int rowptr_[] = {0, 2, 4, 6, 8, 9};
    double values_[] = {2.0, -1.5, 3.0, 1.0, 4.0, -2.0, 1.5, -0.5, 2.5};
    int nnz = rowptr_[n];
    Vector<int> Col(nnz); Vector<long> Row(n+1);
    DVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }

    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];

    A.SetData(n, n, Values, Row, Col);
    // Computation of right hand side.
    DVect x_sol(n), b_vec(n);
    b_vec.Zero();
    x_sol.Fill();
    Mlt(A, x_sol, b_vec);
    x_sol.Zero();

    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    
    success_rowsym_real = CheckSolution(x_sol);
  }

  // We test RowSymSparse (complex numbers).
  {
    Matrix<complex<double>, Symmetric, RowSymSparse> A;
    int colind_[] = {0, 1, 1, 2, 2, 3, 3, 4, 4};
    int rowptr_[] = {0, 2, 4, 6, 8, 9};
    complex<double> values_[] = {cpx(2.0,-4.0), cpx(-1.5,0.5), cpx(3.0,-1.0),
				 cpx(1.0,0.0), cpx(4.0,2.0), cpx(-2.0,1.0),
				 cpx(1.5,-1.0), cpx(-0.5,0.0), cpx(2.5,5.0)};
    int nnz = rowptr_[n];
    Vector<int> Col(nnz); Vector<long> Row(n+1);
    ZVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }

    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];

    A.SetData(n, n, Values, Row, Col);
    // Computation of right hand side
    ZVect x_sol(n), b_vec(n);
    x_sol.Fill(); b_vec.Zero();
    Mlt(A, x_sol, b_vec);
    x_sol.Zero();

    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    success_rowsym_complex = CheckSolution(x_sol);
  }

  // We test RowSparse (real numbers).
  {
    Matrix<double, General, RowSparse> A;
    int colind_[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3, 4, 3, 4};
    int rowptr_[] = {0, 3, 6, 10, 13, 15};
    double values_[] = {2.0, -1.5, 0.4, -1.0, 3.0, 1.0, -2.0,
			1.3, 4.0, -2.0, 0.5, 1.5, -0.5, -0.8, 2.5};

    int nnz = rowptr_[n];
    Vector<int> Col(nnz); Vector<long> Row(n+1);
    DVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }

    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];

    A.SetData(n, n, Values, Row, Col);
    
    // computation of right hand side
    DVect x_sol(n), b_vec(n), b_trans(n);
    b_vec.Zero(); b_trans.Zero();
    x_sol.Fill();
    Mlt(SeldonNoTrans, A, x_sol, b_vec);
    Mlt(SeldonTrans, A, x_sol, b_trans);
    x_sol.Zero();

    // we solve linear system
    Solve(A, x_sol, b_vec, b_trans);
    
    success_row_real = CheckSolution(x_sol);
    successT_row_real = CheckSolution(b_trans);
  }

  // We test RowSparse (complex numbers).
  {
    Matrix<cpx, General, RowSparse> A;
    int colind_[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3, 4, 3, 4};
    int rowptr_[] = {0, 3, 6, 10, 13, 15};
    cpx values_[] = {cpx(2.0,-4.0), cpx(-1.5,0.5), cpx(-0.4,0.2),
		     cpx(1.2,0.2), cpx(3.0,-1.0), cpx(1.0,0.0), cpx(-1.8,0.5),
		     cpx(0.6,0.8), cpx(4.0,2.0), cpx(-2.0,1.0), cpx(1.5,0.4),
		     cpx(1.5,-1.0), cpx(-0.5,0.0), cpx(0.7,-0.5),
		     cpx(2.5,5.0)};

    int nnz = rowptr_[n];
    Vector<int> Col(nnz); Vector<long> Row(n+1);
    ZVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }

    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];

    A.SetData(n, n, Values, Row, Col);
    
    // Computation of right hand side.
    ZVect x_sol(n), b_vec(n), b_trans(n);
    x_sol.Fill();  b_vec.Zero(); b_trans.Zero();
    Mlt(SeldonNoTrans, A, x_sol, b_vec);
    Mlt(SeldonTrans, A, x_sol, b_trans);
    x_sol.Zero();

    // We solve the linear system.
    Solve(A, x_sol, b_vec, b_trans);
    
    success_row_complex = CheckSolution(x_sol);
    successT_row_complex = CheckSolution(b_trans);
  }
  
   // Results for Array...
  bool successA_rowsym_real, successA_rowsym_complex, successA_row_real, successA_row_complex;

  // tranpose
  bool successAT_row_real, successAT_row_complex;

  // We test ArrayRowSymSparse (real numbers).
  {
    // construction of a matrix
    Matrix<double, Symmetric, ArrayRowSymSparse> A(n, n);
    A.Get(0, 0) = 2.0; A.Get(0, 1) = -1.5;
    A.Get(1, 1) = 3.0; A.Get(1, 2) = 1.0;
    A.Get(2, 2) = 4.0; A.Get(2, 3) = -2.0;
    A.Get(3, 3) = 1.5; A.Get(3, 4) = -0.5; A.Get(4, 4) = 2.5;
    
    // Computation of right hand side.
    DVect x_sol(n), b_vec(n); b_vec.Zero();
    x_sol.Fill();
    Mlt(1.0, A, x_sol, b_vec);
    x_sol.Zero();

    // We solve the linear system.
    Solve(A, x_sol, b_vec);

    successA_rowsym_real = CheckSolution(x_sol);
  }

  // We test ArrayRowSymSparse (complex numbers).
  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> A(n, n);
    A.Get(0, 0) = cpx(2.0,-4.0); A.Get(0, 1) = cpx(-1.5,0.5);
    A.Get(1, 1) = cpx(3.0,-1.0); A.Get(1, 2) = cpx(1.0,0.0);
    A.Get(2, 2) = cpx(4.0,2.0); A.Get(2, 3) = cpx(-2.0,1.0);
    A.Get(3, 3) = cpx(1.5,-1.0); A.Get(3, 4) = cpx(-0.5,0.0); A.Get(4, 4) = cpx(2.5,5.0);
    
    // Computation of right hand side
    ZVect x_sol(n), b_vec(n);
    x_sol.Fill(); b_vec.Zero();
    Mlt(A, x_sol, b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    successA_rowsym_complex = CheckSolution(x_sol);
  }
  
  // We test ArrayRowSparse (real numbers).
  {
    Matrix<double, General, ArrayRowSparse> A(n, n);
    A.Get(0, 0) = 2.0; A.Get(0, 1) = -1.5; A.Get(0, 2) = 0.4;
    A.Get(1, 0) = -1.0; A.Get(1, 1) = 3.0; A.Get(1, 2) = 1.0;
    A.Get(2, 0) = -2.0; A.Get(2, 1) = 1.3; A.Get(2, 2) = 4.0; A.Get(2, 3) = -2.0;
    A.Get(3, 2) = 0.5; A.Get(3, 3) = 1.5; A.Get(3, 4) = -0.5;
    A.Get(4, 3) = -0.8; A.Get(4, 4) = 2.5;
    
    // computation of right hand side
    DVect x_sol(n), b_vec(n), b_trans(n);
    x_sol.Fill();  b_vec.Zero(); b_trans.Zero();
    Mlt(SeldonNoTrans, A, x_sol, b_vec);
    Mlt(SeldonTrans, A, x_sol, b_trans);
    x_sol.Zero();

    // we solve linear system
    Solve(A, x_sol, b_vec, b_trans);
    
    successA_row_real = CheckSolution(x_sol);
    successAT_row_real = CheckSolution(b_trans);
  }

  // We test ArrayRowSparse (complex numbers).
  {
    Matrix<cpx, General, ArrayRowSparse> A(n, n);
    A.Get(0, 0) = cpx(2.0,-4.0); A.Get(0, 1) = cpx(-1.5,0.5); A.Get(0, 2) = cpx(-0.4,0.2);
    A.Get(1, 0) = cpx(1.2,0.2); A.Get(1, 1) = cpx(3.0,-1.0); A.Get(1, 2) = cpx(1.0,0.0);
    A.Get(2, 0) = cpx(-1.8,0.5); A.Get(2, 1) = cpx(0.6,0.8); A.Get(2, 2) = cpx(4.0,2.0); A.Get(2, 3) = cpx(-2.0,1.0);
    A.Get(3, 2) = cpx(1.5,0.4); A.Get(3, 3) = cpx(1.5,-1.0); A.Get(3, 4) = cpx(-0.5,0.0);
    A.Get(4, 3) = cpx(0.7,-0.5); A.Get(4, 4) = cpx(2.5,5.0);
    
    // Computation of right hand side.
    ZVect x_sol(n), b_vec(n), b_trans(n);
    x_sol.Fill();  b_vec.Zero(); b_trans.Zero();
    Mlt(SeldonNoTrans, A, x_sol, b_vec);
    Mlt(SeldonTrans, A, x_sol, b_trans);
    x_sol.Zero();

    // We solve the linear system.
    Solve(A, x_sol, b_vec, b_trans);
    
    successA_row_complex = CheckSolution(x_sol);
    successAT_row_complex = CheckSolution(b_trans);
  }

  bool overall_success = true;
  if (!success_rowsym_real)
    {
      cout << "Error during inversion of RowSymSparse real matrix" << endl;
      overall_success = false;
    }

  if (!success_rowsym_complex)
    {
      cout << "Error during inversion of RowSymSparse complex matrix" << endl;
      overall_success = false;
    }

  if (!success_row_real)
    {
      cout << "Error during inversion of RowSparse real matrix" << endl;
      overall_success = false;
    }

  if (!successT_row_real)
    {
      cout << "Error during inversion of transpose RowSparse real matrix" << endl;
      overall_success = false;
    }

  if (!success_row_complex)
    {
      cout << "Error during inversion of RowSparse complex matrix" << endl;
      overall_success = false;
    }

  if (!successT_row_complex)
    {
      cout << "Error during inversion of transpose RowSparse complex matrix" << endl;
      overall_success = false;
    }
  
  if (!successA_rowsym_real)
    {
      cout << "Error during inversion of ArrayRowSymSparse real matrix" << endl;
      overall_success = false;
    }

  if (!successA_rowsym_complex)
    {
      cout << "Error during inversion of ArrayRowSymSparse complex matrix" << endl;
      overall_success = false;
    }

  if (!successA_row_real)
    {
      cout << "Error during inversion of ArrayRowSparse real matrix" << endl;
      overall_success = false;
    }

  if (!successAT_row_real)
    {
      cout << "Error during inversion of transpose ArrayRowSparse real matrix" << endl;
      overall_success = false;
    }

  if (!successA_row_complex)
    {
      cout << "Error during inversion of ArrayRowSparse complex matrix" << endl;
      overall_success = false;
    }

  if (!successAT_row_complex)
    {
      cout << "Error during inversion of transpose ArrayRowSparse complex matrix" << endl;
      overall_success = false;
    }
    
  
  /***************************************************
   * Testing a larger matrix with SparseDirectSolver *
   ***************************************************/
  
  
  {
    Matrix<complex<double>, Symmetric, ArrayRowSymSparse> A;
    Vector<complex<double> > b_vec, x_sol, x_ref;
    
    A.ReadText("test/matrix/MatDisque.dat");
    b_vec.Read("test/matrix/RhsDisque.dat");
    x_ref.Read("test/matrix/SolDisque.dat");
    
    WriteHarwellBoeing(A, "test.rua");
    
    SparseDirectSolver<complex<double> > mat_lu;
    //mat_lu.SelectDirectSolver(mat_lu.SELDON_SOLVER);
    //mat_lu.SelectDirectSolver(mat_lu.MUMPS);
    mat_lu.ShowMessages();
    //mat_lu.SelectOrdering(SparseMatrixOrdering::IDENTITY);
    //mat_lu.SelectOrdering(SparseMatrixOrdering::REVERSE_CUTHILL_MCKEE);
    //mat_lu.SelectOrdering(SparseMatrixOrdering::PORD);
    //mat_lu.SelectOrdering(SparseMatrixOrdering::SCOTCH);
    //mat_lu.SelectOrdering(SparseMatrixOrdering::METIS);
    //mat_lu.SelectOrdering(SparseMatrixOrdering::QAMD);
    mat_lu.Factorize(A);
    
    x_sol = b_vec;
    mat_lu.Solve(x_sol);
    double err;
    bool success = CheckSolution(x_sol, x_ref, err);
    cout << "Error obtained = " << err << endl;
    if (!success)
      {
    	cout << "Error during inversion of ArrayRowSymSparse with SparseDirectSolver" << endl;
    	overall_success = false;
      }
    
    // int test_input; cout << "we wait " << endl; cin >> test_input;
  }
  
  {
    Matrix<complex<double>, General, ArrayRowSparse> A;
    Vector<complex<double> > b_vec, x_sol, x_ref;
    
    A.ReadText("test/matrix/MatDisque.dat");
    b_vec.Read("test/matrix/RhsDisque.dat");
    x_ref.Read("test/matrix/SolDisque.dat");
    
    SparseDirectSolver<complex<double> > mat_lu;
    mat_lu.ShowMessages();
    mat_lu.Factorize(A);
    
    x_sol = b_vec;
    mat_lu.Solve(x_sol);
    double err;
    bool success = CheckSolution(x_sol, x_ref, err);
    cout << "Error obtained = " << err << endl;
    if (!success)
      {
	cout << "Error during inversion of ArrayRowSparse with SparseDirectSolver" << endl;
	overall_success = false;
      }
  }
  
  if (overall_success)
    cout << "All tests successfully completed" << endl;
  else
    return -1;

  return 0;
}
