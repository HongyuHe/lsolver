#ifndef SELDON_FILE_TINY_VECTOR_CXX

#include "TinyVector.hxx"

namespace Seldon
{

  //! intersection of edges [pointA,pointB] and [pt1,pt2]
  /*!
    \param[in] pointA first extremity of the first edge
    \param[in] pointB second extremity of the first edge 
    \param[in] pt1 first extremity of the second edge
    \param[in] pt2 second extremity of the second edge
    \param[out] res intersection
    \param[in] threshold accuracy of points
    \return 0 if no intersection, 1 if intersection and 2 if edges are on the same line
   */
  template<class T>
  int IntersectionEdges(const TinyVector<T,2>& pointA, const TinyVector<T,2>& pointB,
			const TinyVector<T,2>& pt1, const TinyVector<T,2>& pt2,
			TinyVector<T,2>& res, const T& threshold)
  {
    T dx = pointB(0)-pointA(0);
    T dy = pointB(1)-pointA(1);
    T dxj = pt2(0)-pt1(0);
    T dyj = pt2(1)-pt1(1);
    T delta = dx*dyj-dy*dxj;
    T x, y, zero;
    SetComplexZero(zero);
    if (delta != zero)
      {
	x = (dxj*dx*(pointA(1)-pt1(1))+pt1(0)*dyj*dx-pointA(0)*dy*dxj)/delta;
	y = -(dyj*dy*(pointA(0)-pt1(0))+pt1(1)*dxj*dy-pointA(1)*dx*dyj)/delta;
	// DISP(x); DISP(y); DISP((x-pt1(0))*(x-pt2(0)));
	// DISP((y-pt1(1))*(y-pt2(1)));DISP((x-pointA(0))*(x-pointB(0)));
	// DISP((y-pointA(1))*(y-pointB(1)));
	if (((x-pt1(0))*(x-pt2(0)))<=threshold)
	  {
	    if (((y-pt1(1))*(y-pt2(1)))<=threshold)
	      {
		if (((x-pointA(0))*(x-pointB(0)))<=threshold)
		  {
		    if (((y-pointA(1))*(y-pointB(1)))<=threshold)
		      {
			res.Init(x,y);
			return 1;
		      }
		  }
	      }
	  }
      }
    else
      {
	if (abs(dx*(pt1(1)-pointA(1))-dy*(pt1(0)-pointA(0)))<=threshold)
	  {
	    return 2;
	  }
      }  
    return 0;
  }
  
  
  //! intersection of lines [pointA,pointB] and [pt1,pt2]
  /*!
    \param[in] pointA a point of the first line
    \param[in] pointB a point of the first line
    \param[in] pt1 a point of the second line
    \param[in] pt2 a point of the second line
    \param[out] res intersection
    \param[in] threshold accuracy of points
    \return 0 if no intersection, 1 if intersection and 2 if lines are equal
   */
  template<class T>
  int IntersectionDroites(const TinyVector<T,2>& pointA, const TinyVector<T,2>& pointB,
			  const TinyVector<T,2>& pt1, const TinyVector<T,2>& pt2,
			  TinyVector<T,2>& res, const T& threshold)
  {
    T dx = pointB(0)-pointA(0);
    T dy = pointB(1)-pointA(1);
    T dxj = pt2(0)-pt1(0);
    T dyj = pt2(1)-pt1(1);
    T delta = dx*dyj-dy*dxj;
    T x,y, zero;
    SetComplexZero(zero);
    if (abs(delta) > threshold)
      {
	x=(dxj*dx*(pointA(1)-pt1(1))+pt1(0)*dyj*dx-pointA(0)*dy*dxj)/delta;
	if (abs(dx) > abs(dxj))
	  y = dy/dx*(x-pointA(0)) + pointA(1);
	else
	  y = dyj/dxj*(x-pt1(0)) + pt1(1);
	
	res.Init(x,y);
	return 1;
      }
    else
      {
	if (abs(dx*(pt1(1)-pointA(1))-dy*(pt1(0)-pointA(0))) <= threshold)
	  {
	    return 2;
	  }
      }  
    
    return 0;
  }

  
  //! constructing tangent vectors u1 and u2, from normale to a plane
  template<class T>
  void GetVectorPlane(const TinyVector<T,3>& normale,
                      TinyVector<T, 3>& u1, TinyVector<T, 3>& u2)
  {
    Vector<T> coord(3); IVect perm(3);
    perm.Fill();
    for (int i = 0; i < 3; i++)
      coord(i) = abs(normale(i));
    
    Sort(3, coord, perm);
    u1.Zero(); u1(perm(2)) = -normale(perm(1));
    u1(perm(1)) = normale(perm(2));
    
    // le dernier vecteur est obtenu avec un produit vectoriel
    TimesProd(normale, u1, u2);
    
    typename ClassComplexType<T>::Treal one;
    SetComplexOne(one);
    Mlt(one/Norm2(u1), u1); Mlt(one/Norm2(u2), u2);
  }
  
} // end namespace

#define SELDON_FILE_TINY_VECTOR_CXX
#endif
