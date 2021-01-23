#ifndef SELDON_FILE_FUNCTIONS_BASE_CXX

// in this file, we put functions such as Add, Mlt, MltAdd
// with a reduced number of templates in order
// to forbid undesired instantiations of generic functions

/*
  Functions defined in this file:
  
  M X -> Y
  Mlt(M, X, Y)

  alpha M X -> Y
  Mlt(alpha, M, X, Y)

  M X -> Y
  M^T X -> Y
  Mlt(Trans, M, X, Y)

*/

namespace Seldon
{

  /**********************************
   * Functions in Functions_MatVect *
   **********************************/

  
  template<class T0, class Prop0, class Storage0, class Allocator0>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x)
  {
    Vector<T0> xr(x.GetM()), xi(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
	xr(i) = real(x(i));
	xi(i) = imag(x(i));
      }
    
    SolveLuVector(A, pivot, xr);
    SolveLuVector(A, pivot, xi);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T0>(xr(i), xi(i));
  }

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, General, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x)
  {
    Vector<T0> xr(x.GetM()), xi(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
	xr(i) = real(x(i));
	xi(i) = imag(x(i));
      }
    
    SolveLuVector(trans, A, pivot, xr);
    SolveLuVector(trans, A, pivot, xi);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T0>(xr(i), xi(i));
  }

  template<class T0, class Storage0, class Allocator0>
  void SolveLU(const SeldonTranspose& trans,
	       const Matrix<T0, Symmetric, Storage0, Allocator0>& A,
	       const Vector<int>& pivot, Vector<complex<T0> >& x)
  {
    Vector<T0> xr(x.GetM()), xi(x.GetM());
    for (int i = 0; i < x.GetM(); i++)
      {
	xr(i) = real(x(i));
	xi(i) = imag(x(i));
      }
    
    SolveLuVector(A, pivot, xr);
    SolveLuVector(A, pivot, xi);
    
    for (int i = 0; i < x.GetM(); i++)
      x(i) = complex<T0>(xr(i), xi(i));
  }
  
}

#define SELDON_FILE_FUNCTIONS_BASE_CXX
#endif
