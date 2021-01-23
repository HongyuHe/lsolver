#ifndef SELDON_FILE_TINY_MATRIX_HXX

#include "vector/TinyVector.hxx"

#include "TinyMatrixExpression.hxx"

namespace Seldon
{
  //! empty class overloaded for general and symmetric matrices
  template<class T, class Prop, int m_, int n_> 
  class TinyMatrix
  {
  };
  
  template<class Prop, int m, int n, int p>
  class TinyMatrixNode
  {
  };


  template<int p>
  class TinyMatrixLoop;

  template<int p, int q>
  class TinyMatrixDoubleLoop;

  template<int p, int q, int r>
  class TinyMatrixTripleLoop;
  
  /******************
   * General matrix *
   ******************/
  
  //! Class storing small matrices whose number of rows and columns
  //! is known at compilation time
  template <class T, int m, int n>
  class TinyMatrix<T, General, m, n>
    : public TinyMatrixExpression<T, m, n, TinyMatrix<T, General, m, n> >
  {
    template<int p>
    friend class TinyMatrixLoop;
    
    template<int p, int q>
    friend class TinyMatrixDoubleLoop;

    template<int p, int q, int r>
    friend class TinyMatrixTripleLoop;
    
  protected: 
    //! elements of the matrix
    T data_[m*n];
    
  public :
    enum{ size_ = m*n };
    
    TinyMatrix();
    explicit TinyMatrix(int i);
    explicit TinyMatrix(int i, int j);
    
    template<class T0>
    TinyMatrix(const TinyVector<TinyVector<T0, m>, n>& A);

    template<class E>
    TinyMatrix(const TinyMatrixExpression<T, m, n, E>& A);

    // Basic methods.
    static int GetM();
    static int GetN();
    static int GetSize();
    T* GetData();
    
    // product of a matrix by a scalar this = this * a where a is a scalar
    template<class T1>
    TinyMatrix<T, General, m, n> & operator *=(const T1& a );
    
    // sum of matrices this = this + B
    template<class T1, class E>
    TinyMatrix<T, General, m, n> & operator +=(const TinyMatrixExpression<T1, m, n, E>& B);
    
    // sum of matrices this = this - B
    template<class T1, class E>
    TinyMatrix<T, General, m, n> & operator -=(const TinyMatrixExpression<T1, m, n, E>& B);

    // fills matrix with x
    TinyMatrix<T, General, m, n> & operator =(const T& x);

    // *this = expression with matrices
    template<class T0, class E>
    TinyMatrix<T, General, m, n> & operator =(const TinyMatrixExpression<T0, m, n, E>& x);
    
    T& operator()(int i,int j);
    const T& operator()(int i,int j) const;
    
    // Convenient functions
    void Zero();         // Fills the matrix of 0
    void SetIdentity();
    
    template<class T0>
    void SetDiagonal(const T0& a); 
    
    void Fill();
    void FillRand();
    
    template<class T0>
    void Fill(const T0& a); // Fills the matrix of a

    bool IsZero() const;
    
    void Write(const string& file_name) const;
    void Write(ostream& out) const;
    
  };
  
  //! class used to obtain (i,j) from data[p]
  template<int m, int n, int p>
  class TinyMatrixNode<General, m, n, p>
  {
  public :
    enum{ i = p/n, j = p%n };
  };

  
  //! class used to obtain (i,j) from data[p]
  template<int m, int p>
  class TinyMatrixNode<Symmetric, m, m, p>
  {
  public :
    enum{ i = TinyMatrixNode<Symmetric, m, m, p-1>::i
          + (TinyMatrixNode<Symmetric, m, m, p-1>::j==(m-1)),
          j = (TinyMatrixNode<Symmetric, m, m, p-1>::j + 1)
          *(TinyMatrixNode<Symmetric, m, m, p-1>::j < (m-1))
          + (TinyMatrixNode<Symmetric, m, m, p-1>::i + 1)
          *(TinyMatrixNode<Symmetric, m, m, p-1>::j == (m-1))
    };
    
  };
  
  template<int m>
  class TinyMatrixNode<Symmetric, m, m, 0>
  {
  public :
    enum { i = 0, j = 0};
  };

  /*****************************
   * TinyMatrix loop unrolling *
   *****************************/
    
  //! class for triple loop in matrix functions
  template<int p, int q, int r>
  class TinyMatrixTripleLoop
  {
  public :
    template<class T, int m>
    static void AddRow(TinyMatrix<T, General, m, m>& A, const T& val);
    
    template<class T, int m>
    static void ModifyUpperCholesky(TinyMatrix<T, Symmetric, m, m>& A,
                                    T& invVal, T& vloc);

  };

  //! class terminating triple loop in matrix functions
  template<int p, int q>
  class TinyMatrixTripleLoop<p, q, 0>
  {
  public :
    template<class T, int m>
    static inline void AddRow(TinyMatrix<T, General, m, m>& A, const T& val) {}
    
    template<class T, int m>
    static inline void ModifyUpperCholesky(TinyMatrix<T, Symmetric, m, m>& A,
                                           T& invVal, T& vloc) {}

  };
  
  //! class for double loop in matrix functions
  template<int p, int q>
  class TinyMatrixDoubleLoop
  {
  public :
    template<int m, int n, class T0, class T1>
    static void Init(const TinyVector<T0, m>& x, TinyMatrix<T1, General, m, n>& A);

    template<int m, int n, class T, class E>
    static void WriteText(ostream& out, const TinyMatrixExpression<T, m, n, E>& A);
    
    template<int m, int n, class T, class Prop>
    static void Write(ostream& out, const TinyMatrix<T, Prop, m, n>& A);

    template<int m, int n, class T, class E1, class E2>
    static bool IsEqual(const TinyMatrixExpression<T, m, n, E1>& A,
			const TinyMatrixExpression<T, m, n, E2>& B);

    template<int m, int n, class T0, class E0, class T1, class E1, class T2>
    static void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
		    const TinyVectorExpression<T1, n, E1>& x,
		    T2& y);
    
    template<int m, int n, int k, class T0, class E0,
	     class T1, class E1, class T2, class Prop2>
    static void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
		    const TinyMatrixExpression<T1, n, k, E1>& B,
		    TinyMatrix<T2, Prop2, m, k>& C);

    template<int m, int n, class T0, class T1, class E1, class T2, class E2, class T3>
    static void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                            const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A);

    template<int m, int n, class T1, class E1, class T2, class E2, class T3>
    static void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                            const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A);

    template<int m, class T0, class T1, class E1, class T3>
    static void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                            TinyMatrix<T3, Symmetric, m, m>& A);

    template<int m, class T1, class E1, class T3>
    static void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                            TinyMatrix<T3, Symmetric, m, m>& A);
    
    template<class T, int m>
    static void GetMaximumColumn(TinyMatrix<T, General, m, m>& A, int& jmax, T& val);
    
    template<class T, int m>
    static void SwapRow(TinyMatrix<T, General, m, m>& A, int i2, T& val);

    template<class T, int m>
    static void SwapColumn(TinyMatrix<T, General, m, m>& A, int i2, T& val);

    template<class T, int m>
    static void MltRow(TinyMatrix<T, General, m, m>& A, const T& coef);

    template<class T, int m>
    static void PerformElimination(TinyMatrix<T, General, m, m>& A,
				   const T& coef, T& val);
    
    template<class T, int m>
    static void PerformSolve(TinyMatrix<T, General, m, m>& A, T& val);
    
    template<class T, int m>
    static void GetDiagonalCholesky(TinyMatrix<T, Symmetric, m, m>& A, T& val);

    template<class T, int m>
    static void ModifyUpperCholesky(TinyMatrix<T, Symmetric, m, m>& A,
                                    T& invVal, T& vloc);
    
    template<class T, class T2, int m>
    static void SolveCholesky(const class_SeldonNoTrans& trans,
                              const TinyMatrix<T, Symmetric, m, m>& A,
                              TinyVector<T2, m>& x);

    template<class T, class T2, int m>
    static void SolveCholesky(const class_SeldonTrans& trans,
                              const TinyMatrix<T, Symmetric, m, m>& A,
                              TinyVector<T2, m>& x, T2& val);
    
    template<class T, class T2, int m>
    static void MltCholesky(const class_SeldonTrans& trans,
                            const TinyMatrix<T, Symmetric, m, m>& A,
                            TinyVector<T2, m>& x, T2& val);

    template<class T, class T2, int m>
    static void MltCholesky(const class_SeldonNoTrans& trans,
                            const TinyMatrix<T, Symmetric, m, m>& A,
                            TinyVector<T2, m>& x);
    
  };
  
  
  //! class for terminating double loops
  template<int p>
  class TinyMatrixDoubleLoop<p, 0>
  {
  public :
    template<int m, int n, class T0, class T1>
    static inline void Init(const TinyVector<T0, m>& x, TinyMatrix<T1, General, m, n>& A)
    {}

    template<int m, int n, class T, class E>
    static inline void WriteText(ostream& out, const TinyMatrixExpression<T, m, n, E>& A) {}
    
    template<int m, int n, class T, class Prop>
    static inline void Write(ostream& out, const TinyMatrix<T, Prop, m, n>& A) {}

    template<int m, int n, class T, class E1, class E2>
    static inline bool IsEqual(const TinyMatrixExpression<T, m, n, E1>& A,
			       const TinyMatrixExpression<T, m, n, E2>& B) { return true; }

    template<int m, int n, class T0, class E0, class T1, class E1, class T2>
    static inline void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
			   const TinyVectorExpression<T1, n, E1>& x,
			   T2& val) {}
    
    template<int m, int n, int k, class T0, class E0,
	     class T1, class E1, class T2, class Prop2>
    static inline void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
			   const TinyMatrixExpression<T1, n, k, E1>& B,
			   TinyMatrix<T2, Prop2, m, k>& C) {}
    
    template<int m, int n, class T0, class T1, class E1, class T2, class E2, class T3>
    static inline void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                                   const TinyVectorExpression<T2, n, E2>& y,
                                   TinyMatrix<T3, General, m, n>& A)
    {}

    template<int m, int n, class T1, class E1, class T2, class E2, class T3>
    static inline void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                                   const TinyVectorExpression<T2, n, E2>& y,
                                   TinyMatrix<T3, General, m, n>& A)
    {}

    template<int m, class T0, class T1, class E1, class T3>
    static inline void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                                   TinyMatrix<T3, Symmetric, m, m>& A)
    {}

    template<int m, class T1, class E1, class T3>
    static inline void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                                   TinyMatrix<T3, Symmetric, m, m>& A)
    {}

    template<class T, int m>
    static inline void GetMaximumColumn(TinyMatrix<T, General, m, m>& A, int& jmax, T& val) {}
    
    template<class T, int m>
    static inline void SwapRow(TinyMatrix<T, General, m, m>& A, int i2, T& val) {}

    template<class T, int m>
    static inline void SwapColumn(TinyMatrix<T, General, m, m>& A, int i2, T& val) {}

    template<class T, int m>
    static inline void MltRow(TinyMatrix<T, General, m, m>& A, const T& coef) {}

    template<class T, int m>
    static inline void PerformElimination(TinyMatrix<T, General, m, m>& A,
					  const T& coef, T& val) {}

    template<class T, int m>
    static inline void PerformSolve(TinyMatrix<T, General, m, m>& A, T& val) {}
    
    template<class T, int m>
    static inline void GetDiagonalCholesky(TinyMatrix<T, Symmetric, m, m>& A, T& val) {}
    
    template<class T, int m>
    static inline void ModifyUpperCholesky(TinyMatrix<T, Symmetric, m, m>& A,
					   T& invVal, T& vloc) {}
    
    template<class T, class T2, int m>
    static inline void SolveCholesky(const class_SeldonNoTrans& trans,
				     const TinyMatrix<T, Symmetric, m, m>& A,
				     TinyVector<T2, m>& x) {}

    template<class T, class T2, int m>
    static inline void SolveCholesky(const class_SeldonTrans& trans,
                                     const TinyMatrix<T, Symmetric, m, m>& A,
                                     TinyVector<T2, m>& x, T2& val) {}
    
    template<class T, class T2, int m>
    static inline void MltCholesky(const class_SeldonTrans& trans,
                                   const TinyMatrix<T, Symmetric, m, m>& A,
                                   TinyVector<T2, m>& x, T2& val) {}

    template<class T, class T2, int m>
    static inline void MltCholesky(const class_SeldonNoTrans& trans,
                                   const TinyMatrix<T, Symmetric, m, m>& A,
                                   TinyVector<T2, m>& x) {}

  };
  
  //! class for simple loops for matrix operations
  template<int p>
  class TinyMatrixLoop
  {
  public :
    template<int m, int n, class T0, class Prop>
    static void Zero(TinyMatrix<T0, Prop, m, n>& A);
    
    template<int m, int n, class T0, class T1>
    static void Init(const TinyVector<TinyVector<T0, m>, n>& x,
                     TinyMatrix<T1, General, m, n>& A);

    template<int m, int n, class T0, class Prop, class T1>
    static void MltScal(const T0& alpha, TinyMatrix<T1, Prop, m, n>& A);
        
    template<int m, int n, class T1, class E, class T0, class Prop>
    static void Copy(const TinyMatrixExpression<T1, m, n, E>& x,
		     TinyMatrix<T0, Prop, m, n>& y);

    template<int m, int n, class T1, class E, class T0, class Prop>
    static void AddCopy(const TinyMatrixExpression<T1, m, n, E>& x,
			TinyMatrix<T0, Prop, m, n>& y);

    template<int m, int n, class T1, class E, class T0, class Prop>
    static void DiffCopy(const TinyMatrixExpression<T1, m, n, E>& x,
			 TinyMatrix<T0, Prop, m, n>& y);
    
    template<int m, int n, class T, class Prop>
    static void SetIdentity(TinyMatrix<T, Prop, m, n>& A);

    template<int m, int n, class T, class Prop, class T0>
    static void SetDiagonal(TinyMatrix<T, Prop, m, n>& A, const T0& diag);
    
    template<int m, int n, class T, class Prop>
    static void Fill(TinyMatrix<T, Prop, m, n>& A);

    template<int m, int n, class T, class Prop>
    static void FillRand(TinyMatrix<T, Prop, m, n>& A);
    
    template<int m, int n, class T, class Prop, class T0>
    static void Fill(TinyMatrix<T, Prop, m, n>& A, const T0& alpha);

    template<int m, int n, class Prop, class T>
    static bool IsZero(const TinyMatrix<T, Prop, m, n>& A);
        
    template<int m, int n, class T, class E>
    static void WriteText(ostream& out, const TinyMatrixExpression<T, m, n, E>& A);

    template<int m, int n, class T, class Prop>
    static void Write(ostream& out, const TinyMatrix<T, Prop, m, n>& A);

    template<int m, int n, class T, class E1, class E2>
    static bool IsEqual(const TinyMatrixExpression<T, m, n, E1>& A,
                        const TinyMatrixExpression<T, m, n, E2>& B);
    
    template<int m, int n, class T0, class E0, class T1, class E1, class T2>
    static void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
		    const TinyVectorExpression<T1, n, E1>& x,
		    TinyVector<T2, m>& y);
    
    template<int m, int n, class T0, class E0, class T1, class E1, class T2>
    static void MltAdd(const TinyMatrixExpression<T0, m, n, E0>& A,
                       const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y);
    
    template<int m, int n, class T0, class E0, class T1, class E1, class T2, class T3>
    static void MltAdd(const T3& alpha, const TinyMatrixExpression<T0, m, n, E0>& A,
		       const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y);
    
    template<int m, int n, int k, class T0, class E0,
	     class T1, class E1, class T2, class Prop2>
    static void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
		    const TinyMatrixExpression<T1, n, k, E1>& B,
		    TinyMatrix<T2, Prop2, m, k>& C);    
    
    template<int m, int n, class T0, class T1, class E1, class T2, class E2, class T3>
    static void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                            const TinyVectorExpression<T2, n, E2>& y,
			    TinyMatrix<T3, General, m, n>& A);

    template<int m, int n, class T1, class E1, class T2, class E2, class T3>
    static void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                            const TinyVectorExpression<T2, n, E2>& y,
			    TinyMatrix<T3, General, m, n>& A);

    template<int m, class T0, class T1, class E1, class T3>
    static void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                            TinyMatrix<T3, Symmetric, m, m>& A);

    template<int m, class T1, class E1, class T3>
    static void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                            TinyMatrix<T3, Symmetric, m, m>& A);
    
    template<int m, int n, class T1, class E1>
    static void GetCol(const TinyMatrixExpression<T1, m, n, E1>& A,
                       int k, TinyVector<T1, m>& x);
    
    template<int m, int n, class T0, class E0, class T1>
    static void GetRow(const TinyMatrixExpression<T0, m, n, E0>& A,
                       int k, TinyVector<T1, n>& x);
    
    template<int m, int n, class T1, class E1, class Prop>
    static void SetCol(const TinyVectorExpression<T1, m, E1>& x,
                       int k, TinyMatrix<T1, Prop, m, n>& A);
    
    template<int m, int n, class T1, class E1, class Prop>
    static void SetRow(const TinyVectorExpression<T1, n, E1>& x,
                       int k, TinyMatrix<T1, Prop, m, n>& A);
    
    template<int m, int n, class T, class Prop, class T0>
    static void GetMaxAbs(const TinyMatrix<T, Prop, m, n>& A, T0& amax);
        
    template<class T, int m>
    static void PivotGauss(TinyMatrix<T, General, m, m>& A, TinyVector<int, m>& pivot);
    
    template<class T, int m>
    static void SolveUpper(TinyMatrix<T, General, m, m>& A);    

    template<class T, int m>
    static void PermuteColumn(TinyMatrix<T, General, m, m>& A, const TinyVector<int, m>& pivot);

    template<class T, int m>
    static void GetCholesky(TinyMatrix<T, Symmetric, m, m>& A);
    
    template<class T, class T2, int m>
    static void SolveCholesky(const class_SeldonNoTrans& trans,
                              const TinyMatrix<T, Symmetric, m, m>& A,
                              TinyVector<T2, m>& x);

    template<class T, class T2, int m>
    static void SolveCholesky(const class_SeldonTrans& trans,
                              const TinyMatrix<T, Symmetric, m, m>& A,
                              TinyVector<T2, m>& x);

    template<class T, class T2, int m>
    static void MltCholesky(const class_SeldonTrans& trans,
                            const TinyMatrix<T, Symmetric, m, m>& A,
                            TinyVector<T2, m>& x);
    
    template<class T, class T2, int m>
    static void MltCholesky(const class_SeldonNoTrans& trans,
                            const TinyMatrix<T, Symmetric, m, m>& A,
                            TinyVector<T2, m>& x);

  };
  
  //! class terminating simple loops for matrix operations
  template<>
  class TinyMatrixLoop<0>
  {
  public :
    template<int m, int n, class T0, class Prop>
    static inline void Zero(TinyMatrix<T0, Prop, m, n>& A) {}
    
    template<int m, int n, class T0, class T1>
    static inline void Init(const TinyVector<TinyVector<T0, m>, n>& x,
                            TinyMatrix<T1, General, m, n>& A) {}

    template<int m, int n, class T0, class Prop, class T1>
    static inline void MltScal(const T0& alpha, TinyMatrix<T1, Prop, m, n>& A) {}
    
    template<int m, int n, class T1, class E, class T0, class Prop>
    static inline void Copy(const TinyMatrixExpression<T1, m, n, E>& x,
			    TinyMatrix<T0, Prop, m, n>& y) {}

    template<int m, int n, class T1, class E, class T0, class Prop>
    static inline void AddCopy(const TinyMatrixExpression<T1, m, n, E>& x,
			       TinyMatrix<T0, Prop, m, n>& y) {}

    template<int m, int n, class T1, class E, class T0, class Prop>
    static inline void DiffCopy(const TinyMatrixExpression<T1, m, n, E>& x,
				TinyMatrix<T0, Prop, m, n>& y) {}

    template<int m, int n, class T, class Prop>
    static inline void SetIdentity(TinyMatrix<T, Prop, m, n>& A) {}

    template<int m, int n, class T, class Prop, class T0>
    static inline void SetDiagonal(TinyMatrix<T, Prop, m, n>& A, const T0& alpha) {}
    
    template<int m, int n, class T, class Prop>
    static inline void Fill(TinyMatrix<T, Prop, m, n>& A) {}

    template<int m, int n, class T, class Prop>
    static inline void FillRand(TinyMatrix<T, Prop, m, n>& A) {}
    
    template<int m, int n, class T, class Prop, class T0>
    static inline void Fill(TinyMatrix<T, Prop, m, n>& A, const T0& alpha) {}

    template<int m, int n, class Prop, class T>
    static inline bool IsZero(const TinyMatrix<T, Prop, m, n>& A) { return true; }
    
    template<int m, int n, class T, class E>
    static inline void WriteText(ostream& out, const TinyMatrixExpression<T, m, n, E>& A) {}
    
    template<int m, int n, class T, class Prop>
    static inline void Write(ostream& out, const TinyMatrix<T, Prop, m, n>& A) {}

    template<int m, int n, class T, class E1, class E2>
    static inline bool IsEqual(const TinyMatrixExpression<T, m, n, E1>& A,
                               const TinyMatrixExpression<T, m, n, E2>& B)
    { return true; }
    
    template<int m, int n, class T0, class E0, class T1, class E1, class T2>
    static inline void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
			   const TinyVectorExpression<T1, n, E1>& x,
			   TinyVector<T2, m>& y) {}
    
    template<int m, int n, class T0, class E0, class T1, class E1, class T2>
    static inline void MltAdd(const TinyMatrixExpression<T0, m, n, E0>& A,
			      const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y) {}
    
    template<int m, int n, class T0, class E0, class T1, class E1, class T2, class T3>
    static inline void MltAdd(const T3& alpha, const TinyMatrixExpression<T0, m, n, E0>& A,
			      const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y) {}
    
    template<int m, int n, int k, class T0, class E0,
	     class T1, class E1, class T2, class Prop2>
    static inline void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
			   const TinyMatrixExpression<T1, n, k, E1>& B,
			   TinyMatrix<T2, Prop2, m, k>& C) {}

    template<int m, int n, class T0, class T1, class E1, class T2, class E2, class T3>
    static inline void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                                   const TinyVectorExpression<T2, n, E2>& y,
                                   TinyMatrix<T3, General, m, n>& A) {}

    template<int m, int n, class T1, class E1, class T2, class E2, class T3>
    static inline void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                                   const TinyVectorExpression<T2, n, E2>& y,
                                   TinyMatrix<T3, General, m, n>& A) {}

    template<int m, class T0, class T1, class E1, class T3>
    static inline void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                                   TinyMatrix<T3, Symmetric, m, m>& A) {}

    template<int m, class T1, class E1, class T3>
    static inline void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                                   TinyMatrix<T3, Symmetric, m, m>& A) {}
    
    template<int m, int n, class T1, class E1>
    static inline void GetCol(const TinyMatrixExpression<T1, m, n, E1>& A,
                              int k, TinyVector<T1, m>& x ) {}
    
    template<int m, int n, class T0, class E0, class T1>
    static inline void GetRow(const TinyMatrixExpression<T0, m, n, E0>& A,
                              int k, TinyVector<T1, n>& x ) {}
    
    template<int m, int n, class T1, class E1, class Prop>
    static inline void SetCol(const TinyVectorExpression<T1, m, E1>& x ,
                              int k,  TinyMatrix<T1, Prop, m, n>& A) {}
    
    template<int m, int n, class T1, class E1, class Prop>
    static inline void SetRow(const TinyVectorExpression<T1, n, E1>& x ,
                              int k,  TinyMatrix<T1, Prop, m, n>& A) {}

    template<int m, int n, class T, class Prop, class T0>
    static inline void GetMaxAbs(const TinyMatrix<T, Prop, m, n>& A, T0& amax)
    {}
    
    template<class T, int m>
    static inline void PivotGauss(TinyMatrix<T, General, m, m>& A,
                                  TinyVector<int, m>& pivot) {}

    template<class T, int m>
    static inline void SolveUpper(TinyMatrix<T, General, m, m>& A) {}
    
    template<class T, int m>
    static inline void PermuteColumn(TinyMatrix<T, General, m, m>& A,
                                     const TinyVector<int, m>& pivot) {}
    
    template<class T, int m>
    static inline void GetCholesky(TinyMatrix<T, Symmetric, m, m>& A) {}

    template<class T, class T2, int m>
    static inline void SolveCholesky(const class_SeldonNoTrans& trans,
                              const TinyMatrix<T, Symmetric, m, m>& A,
                              TinyVector<T2, m>& x) {}

    template<class T, class T2, int m>
    static inline void SolveCholesky(const class_SeldonTrans& trans,
                              const TinyMatrix<T, Symmetric, m, m>& A,
                              TinyVector<T2, m>& x) {}
    
    template<class T, class T2, int m>
    static inline void MltCholesky(const class_SeldonTrans& trans,
                                   const TinyMatrix<T, Symmetric, m, m>& A,
                                   TinyVector<T2, m>& x) {}
    
    template<class T, class T2, int m>
    static inline void MltCholesky(const class_SeldonNoTrans& trans,
                                   const TinyMatrix<T, Symmetric, m, m>& A,
                                   TinyVector<T2, m>& x) {}

  };
  

  /********************
   * Symmetric matrix *
   ********************/
    
  //! class storing tiny small matrices
  template <class T, int m>
  class TinyMatrix<T, Symmetric, m, m>
    : public TinyMatrixExpression<T, m, m, TinyMatrix<T, Symmetric, m, m> >
  {
    template<int p>
    friend class TinyMatrixLoop;

    template<int p, int q>
    friend class TinyMatrixDoubleLoop;

    template<int p, int q, int r>
    class TinyMatrixTripleLoop;
    
  protected :
    //! elements stored in the matrix
    T data_[m*(m+1)/2];
    
  public :
    enum{ size_ = m*(m+1)/2 };
    
    TinyMatrix();
    explicit TinyMatrix(int i);
    explicit TinyMatrix(int i, int j);

    template<class E>
    TinyMatrix(const TinyMatrixExpression<T, m, m, E>& A);
    
    // Basic methods.
    static int GetM();
    static int GetN();
    static int GetSize();
    T* GetData();
    
    // product of a matrix by a scalar this = this * a where a is a scalar
    template<class T1>
    TinyMatrix<T, Symmetric, m, m> & operator *=(const T1& a );
    
    // sum of matrices res = this + B where B is a symmetric matrix
    template<class T1, class E>
    TinyMatrix<T, Symmetric, m, m> & operator +=(const TinyMatrixExpression<T1, m, m, E>& B);
    
    // sum of matrices res = this - B where B is a symmetric matrix
    template<class T1, class E>
    TinyMatrix<T, Symmetric, m, m> & operator -=(const TinyMatrixExpression<T1, m, m, E>& B);

    // fills matrix with x
    TinyMatrix<T, Symmetric, m, m> & operator =(const T& x);

    // this = expresion with matrices
    template<class T0, class E>
    TinyMatrix<T, Symmetric, m, m> & operator =(const TinyMatrixExpression<T0, m, m, E>& x);

    T& operator()(int i,int j);
    const T& operator()(int i,int j) const;

    // Convenient functions
    void Zero();         // Fills the matrix of 0
    void SetIdentity();
    
    template<class T0>
    void SetDiagonal(const T0& alpha);
    
    void Fill();
    void FillRand();
    
    template<class T1>
    void Fill(const T1& a); // Fills the matrix of a

    bool IsZero() const;
    
    void WriteText(const string& file_name) const;
    void Write(const string& file_name) const;
    void Write(ostream& out) const;
    
  };  
  
  
  /*************
   * Operators *
   *************/

  template<class T, int m, int n, class E1, class E2>
  bool operator==(const TinyMatrixExpression<T, m, n, E1> & u,
		  const TinyMatrixExpression<T, m, n, E2> & v);
  
  template<class T, int m, int n, class E1, class E2>
  bool operator!=(const TinyMatrixExpression<T, m, n, E1> & u,
		  const TinyMatrixExpression<T, m, n, E2> & v);    

  template <class T, int m, int n, class E>
  ostream& operator <<(ostream& out, const TinyMatrixExpression<T, m, n, E> & A);
  

  /****************************
   * Matrix-vector operations *
   ****************************/

  // returns A*x
  template<class T, int m, int n, class E, class T1, class E1> TinyVector<T1, m>
  dot(const TinyMatrixExpression<T, m, n, E>& A,
      const TinyVectorExpression<T1, n, E1>& x);
  
  // returns A*B
  template<class T, int m, int n, class E, class T1, int k, class E1>
  TinyMatrix<T, General, m, k>
  dot(const TinyMatrixExpression<T, m, n, E>& A,
      const TinyMatrixExpression<T1, n, k, E1> & B); 
  
  // product of a matrix by a vector y = A * x
  template<class T0, class E0, class T1, class E1, class T2, int m, int n>
  void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
	   const TinyVectorExpression<T1, n, E1>& x, TinyVector<T2, m>& y);
  
  // product of a matrix by a vector y = y + A*x
  template<class T1, class E1, class T2, class E2, class T3, int m, int n>
  void MltAdd(const TinyMatrixExpression<T1, m, n, E1>& A,
              const TinyVectorExpression<T2, n, E2>& x, TinyVector<T3, m>& y);
  
  // product of a matrix by a vector y = y + alpha*A*x
  template<class T0, class T1, class E1, class T2, class E2,
	   class T3, int m, int n>
  void MltAdd(const T0& alpha, const TinyMatrixExpression<T1, m, n, E1>& A,
              const TinyVectorExpression<T2, n, E2>& x, TinyVector<T3, m>& y);

  template<class T0, class E0, class T1, class E1, class T2, int m, int n>
  void MltTrans(const TinyMatrixExpression<T0, m, n, E0>& A,
		const TinyVectorExpression<T1, m, E1>& x, TinyVector<T2, n>& y);

  template<class T0, class E0, class T1, class E1, class T2, int m, int n>
  void Mlt(const class_SeldonTrans&,
	   const TinyMatrixExpression<T0, m, n, E0>& A,
	   const TinyVectorExpression<T1, m, E1>& x, TinyVector<T2, n>& y);
  
  template<int m, int n, class T0, class T1, class E1,
	   class T2, class E2, class T3>
  void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
                   const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A);

  template<int m, int n, class T1, class E1, class T2, class E2, class T3>
  void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
                   const TinyVectorExpression<T2, n, E2>& y, TinyMatrix<T3, General, m, n>& A);

  template<int m, class T0, class T1, class E1, class T3>
  void Rank1Update(const T0& alpha, const TinyVectorExpression<T1, m, E1>& x,
		   TinyMatrix<T3, Symmetric, m, m>& A);

  template<int m, class T1, class E1, class T3>
  void Rank1Matrix(const TinyVectorExpression<T1, m, E1>& x,
		   TinyMatrix<T3, Symmetric, m, m>& A);
  
  template<int m, int n, class T1, class E1>
  void GetCol(const TinyMatrixExpression<T1, m, n, E1>& A, int k, TinyVector<T1, m>& x);
    
  template<int m, int n, class T0, class E0, class T1>
  void GetRow(const TinyMatrixExpression<T0, m, n, E0>& A, int k, TinyVector<T1, n>& x);
  
  template<int m, int n, class T1, class E1, class Prop>
  void SetCol(const TinyVectorExpression<T1, m, E1>& x, int k, TinyMatrix<T1, Prop, m, n>& A);
  
  template<int m, int n, class T1, class E1, class Prop>
  void SetRow(const TinyVectorExpression<T1, n, E1>& x, int k, TinyMatrix<T1, Prop, m, n>& A);
  
  
  /****************************
   * Matrix-Matrix operations *
   ****************************/
  
  // B = A
  template<class T, int m, int n, class E, class Prop>
  void Copy(const TinyMatrixExpression<T, m, n, E>& A, TinyMatrix<T, Prop, m, n>& B);
  
  // C = A + B
  template<class T0, class E0, class T1, class E1,
	   class T2, class Prop2, int m, int n>
  void Add(const TinyMatrixExpression<T0, m, n, E0>& A,
           const TinyMatrixExpression<T1, m, n, E1>& B, TinyMatrix<T2, Prop2, m, n>& C);
  
  // B = B + alpha*A 
  template<class T1, class E1, class T2, class Prop2, int m, int n>
  void Add(const T1& alpha, const TinyMatrixExpression<T1, m, n, E1>& A,
           TinyMatrix<T2, Prop2, m, n>& B);
  
  // B = B + A 
  template<class T0, class E0, class T1, class Prop1, int m, int n>
  void Add(const TinyMatrixExpression<T0, m, n, E0>& A, TinyMatrix<T1, Prop1, m, n>& B);
  
  // product of a matrix by a scalar A = alpha * A
  template<class T0, class T1, class Prop, int m, int n>
  void Mlt(const T0& alpha, TinyMatrix<T1, General, m, n>& A);
  
  // product of a matrix by a matrix C = A * B
  template<class T0, class E0, class T1, class E1,
	   class T2, class Prop2, int m, int n, int k>
  void Mlt(const TinyMatrixExpression<T0, m, n, E0>& A,
	   const TinyMatrixExpression<T1, n, k, E1>& B,
	   TinyMatrix<T2, Prop2, m, k>& C);

  // product of a matrix by a matrix C = A * B^t
  template<class T0, class E0, class T1, class E1,
	   class T2, class Prop2, int m, int n, int k>
  void MltTrans(const TinyMatrixExpression<T0, m, n, E0>& A,
                const TinyMatrixExpression<T1, k, n, E1>& B, TinyMatrix<T2, Prop2, m, k>& C);
  
  template<class T3, class T0, class E0, class T1, class E1,
	   class T2, class Prop2, class T4, int m, int n, int k>
  void MltAdd(const T3& alpha, const class_SeldonTrans&,
	      const TinyMatrixExpression<T0, n, m, E0>& A,
	      const class_SeldonNoTrans&, const TinyMatrixExpression<T1, n, k, E1>& B,
	      const T4& beta, TinyMatrix<T2, Prop2, m, k>& C);
  
  template<class T, int m, int n, class E>
  void Transpose(const TinyMatrixExpression<T, m, n, E> & A, TinyMatrix<T, General, n, m> & B);
  
  template<class T, int m>
  void Transpose(TinyMatrix<T, General, m, m> & B);

  template<class T, class Prop, int m, int n>
  typename ClassComplexType<T>::Treal MaxAbs(const TinyMatrix<T, Prop, m, n>& A);

  /***************************************
   * Specific functions for 1x1 matrices *
   ***************************************/


  template<class T, class Prop>
  T Det(const TinyMatrix<T, Prop, 1, 1> & A);

  template<class T>
  void GetInverse(const TinyMatrix<T, General, 1, 1> & A,
		  TinyMatrix<T, General, 1, 1> & B);

  template<class T>
  void GetInverse(const TinyMatrix<T, Symmetric, 1, 1> & A,
		  TinyMatrix<T, Symmetric, 1, 1> & B);

  template<class T>
  void GetInverse(TinyMatrix<T, General, 1, 1> & B);

  template<class T>
  void GetInverse(TinyMatrix<T, Symmetric, 1, 1> & B);

  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 1>& n, TinyMatrix<T1, Symmetric, 1, 1>& P);

  template<class T0, class T1>
  void GetNormalProjector(const T0& n, TinyMatrix<T1, Symmetric, 1, 1>& P);

  
  /***************************************
   * Specific functions for 2x2 matrices *
   ***************************************/
  
  
  // determinant of a matrix 2x2
  template<class T, class Prop>
  T Det(const TinyMatrix<T, Prop, 2, 2> & A);

  template<class T0, class T1>
  void GetTangentialProjector(const TinyVector<T0, 2>& n,
			      TinyMatrix<T1, Symmetric, 2, 2>& P);
  
  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 2>& n, TinyMatrix<T1, Symmetric, 2, 2>& P);
  
  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 2>& n, TinyMatrix<T1, General, 2, 2>& P);
  
  // non-inline functions
  template<class T>
  void GetInverse(const TinyMatrix<T, General, 2, 2> & A, TinyMatrix<T, General, 2, 2> & B);
  
  template<class T>
  void GetInverse(TinyMatrix<T, General, 2, 2> & B);
  
  template<class T>
  void GetInverse(const TinyMatrix<T, Symmetric, 2, 2> & A,
                  TinyMatrix<T, Symmetric, 2, 2> & B);
  
  template<class T>
  void GetInverse(TinyMatrix<T, Symmetric, 2, 2> & B);
  

  template<class T>
  void GetEigenvalues(TinyMatrix<T, General, 2, 2>& A,
		      TinyVector<T, 2> & LambdaR, TinyVector<T, 2>& LambdaI);
  
  template<class T>
  void GetEigenvalues(TinyMatrix<complex<T>, General, 2, 2>& A,
		      TinyVector<complex<T>, 2> & Lambda);
  
  template<class T>
  void GetSquareRoot(TinyMatrix<T, Symmetric, 2, 2>& A);
  
  
  /***************************************
   * Specific functions for 3x3 matrices *
   ***************************************/
  

  // non-inline functions
  template<class T, class Prop>
  T Det(const TinyMatrix<T, Prop, 3, 3> & A);

  template<class T>
  void GetInverse(const TinyMatrix<T, General, 3, 3> & A,
                  TinyMatrix<T, General, 3, 3> & B);
  
  template<class T>
  void GetInverse(TinyMatrix<T, General, 3, 3> & B);
  
  template<class T>
  void GetInverse(const TinyMatrix<T, Symmetric, 3, 3> & A,
                  TinyMatrix<T, Symmetric, 3, 3> & B);
  
  template<class T>
  void GetInverse(TinyMatrix<T, Symmetric, 3, 3> & B);

  template<class T0, class T1>
  void GetTangentialProjector(const TinyVector<T0, 3>& n,
			      TinyMatrix<T1, Symmetric, 3, 3>& P);
  
  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 3>& n,
			  TinyMatrix<T1, Symmetric, 3, 3>& P);

  template<class T0, class T1>
  void GetNormalProjector(const TinyVector<T0, 3>& n,
			  TinyMatrix<T1, General, 3, 3>& P);

    
  /*******************
   * Other functions *
   *******************/

  
  template<class T0, class T1, int m>
  void GetTangentialProjector(const TinyVector<T0, m>& n,
			      TinyMatrix<T1, Symmetric, m, m>& P);  

  template<class T0, class T1, int m>
  void GetNormalProjector(const TinyVector<T0, m>& n,
			  TinyMatrix<T1, Symmetric, m, m>& P);

  template<class T, class Prop, int m, int n>
  void FillZero(TinyMatrix<T, Prop, m, n>& X);
  
  //! res = A . B, where the scalar product is performed between columns of A and B
  template<class T, int m, int n>
  void DotProdCol(const TinyMatrix<T, General, m, n>& A,
		  const TinyMatrix<T, General, m, n>& B, TinyVector<T, n>& res);
  
  // non-inline functions
  template<class T, class Property, int m, int n>
  typename ClassComplexType<T>::Treal
  Norm2_Column(const TinyMatrix<T, Property, m, n>& A,
	       int first_row, int index_col);

  template<class T, class Property, class Storage, class Allocator>
  typename ClassComplexType<T>::Treal
  Norm2_Column(const Matrix<T, Property, Storage, Allocator>& A,
	       int first_row, int index_col);
  
  template<class T, int m>
  void GetInverse(TinyMatrix<T, General, m, m>& A);
  
  template<class T, int m>
  void GetInverse(const TinyMatrix<T, General, m, m>& A,
		  TinyMatrix<T, General, m, m>& B);
  
  template<class T, int m>
  void GetInverse(TinyMatrix<T, Symmetric, m, m> & A);

  template<class T, int m>
  void GetInverse(const TinyMatrix<T, Symmetric, m, m>& A,
		  TinyMatrix<T, Symmetric, m, m>& B);

  template<class T, int m>
  void GetCholesky(TinyMatrix<T, Symmetric, m, m>& A);
  
  template<class T, class T2, int m>
  void SolveCholesky(const class_SeldonNoTrans& trans,
                     const TinyMatrix<T, Symmetric, m, m>& A,
                     TinyVector<T2, m>& x);
  
  template<class T, class T2, int m>
  void SolveCholesky(const class_SeldonTrans& trans,
                     const TinyMatrix<T, Symmetric, m, m>& A,
                     TinyVector<T2, m>& x);
  
  template<class T, class T2, int m>
  void MltCholesky(const class_SeldonNoTrans& trans,
                   const TinyMatrix<T, Symmetric, m, m>& A,
                   TinyVector<T2, m>& x);
  
  template<class T, class T2, int m>
  void MltCholesky(const class_SeldonTrans& trans,
                     const TinyMatrix<T, Symmetric, m, m>& A,
                   TinyVector<T2, m>& x);
  
  template<int m, class T>
  void GetEigenvaluesEigenvectors(TinyMatrix<T, Symmetric, m, m>& A,
				  TinyVector<T, m>& w,
				  TinyMatrix<T, General, m, m>& z);
  
  template<int m, class T>
  void GetEigenvalues(TinyMatrix<T, Symmetric, m, m>& A,
		      TinyVector<T, m>& w);

  template<class T, int m>
  void GetSquareRoot(TinyMatrix<T, Symmetric, m, m>& A);
    
} // end namespace


#define SELDON_FILE_TINY_MATRIX_HXX
#endif

