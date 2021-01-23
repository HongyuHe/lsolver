#ifndef SELDON_FILE_MY_MULTIVEC_HPP
#define SELDON_FILE_MY_MULTIVEC_HPP

namespace Anasazi
{

  //! Simple example of a user's defined Anasazi::MultiVec class.
  /*! 
   * This is a simple, single processor example of user's defined
   * MultiVec-derived class. The class is templated with ScalarType;
   * possible choices are, for example, "float", "double", or
   * "complex<double>".
   *
   * \author Oscar Chinallato (ETHZ/ICOS) and Marzio Sala (ETHZ/COLAB) 
   *
   * \date Last modified on 01-Nov-05
   */
  template <class ScalarType>
  class MyMultiVec : public Anasazi::MultiVec<ScalarType>
  {
  public:
    

    MyMultiVec(const int Length, const int NumberVecs);
    MyMultiVec(const int Length, const std::vector<ScalarType*>& rhs);
    MyMultiVec(const MyMultiVec& rhs);

    ~MyMultiVec();

    MyMultiVec<ScalarType>* Clone(const int NumberVecs) const;
    MyMultiVec<ScalarType>* CloneCopy() const;
    MyMultiVec<ScalarType>* CloneCopy(const std::vector< int > &index) const;
    MyMultiVec<ScalarType>* CloneViewNonConst(const std::vector< int > &index);
    const MyMultiVec<ScalarType>* CloneView(const std::vector< int > &index) const;
    
    int GetVecLength () const;
    int GetNumberVecs () const;

    void MvTimesMatAddMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType> &A, 
			  const Teuchos::SerialDenseMatrix<int, ScalarType> &B, 
			  ScalarType beta);

    void MvAddMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
		  ScalarType beta,  const Anasazi::MultiVec<ScalarType>& B);

    void MvTransMv (ScalarType alpha, const Anasazi::MultiVec<ScalarType>& A, 
		    Teuchos::SerialDenseMatrix< int, ScalarType >& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		    , Anasazi::ConjType conj
#endif
		    ) const;

    void MvDot (const Anasazi::MultiVec<ScalarType>& A, std::vector<ScalarType> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
		, Anasazi::ConjType conj
#endif
		) const;
    
    void MvNorm ( std::vector<typename Teuchos::
		  ScalarTraits<ScalarType>::magnitudeType> &normvec ) const;

    void SetBlock (const Anasazi::MultiVec<ScalarType>& A, 
		   const std::vector<int> &index);

    void MvScale( ScalarType alpha );

    void MvScale( const std::vector<ScalarType>& alpha );

    void  MvRandom ();
    void  MvInit (ScalarType alpha);
    void MvPrint (std::ostream &os) const;
    
    ScalarType& operator()(const int i, const int j);
    const ScalarType& operator()(const int i, const int j) const;
    ScalarType* operator[](int v);
    ScalarType* operator[](int v) const;
    
  private:
    void Check();
    
    //! Length of the vectors
    const int Length_;
    //! Number of multi-vectors
    const int NumberVecs_;
    //! Pointers to the storage of the vectors.
    std::vector<ScalarType*> data_;
    //! If \c true, then this object owns the vectors and must free them in dtor.
    std::vector<bool> ownership_;
    
  };
  
  template <typename ScalarType>
  std::ostream& operator<<(std::ostream& os, const MyMultiVec<ScalarType>& Obj);
  
  
}

#endif // MY_MULTIVECTOR_HPP
