// -*-c++-*-
/// \file kernel.h
/// \brief Interface to kernel classes
/// \author Artem Sokolov

#ifndef KERNEL_H__INCLUDED
#define KERNEL_H__INCLUDED

#include "types.h"
#include "params.h"

#include <vector>
#include <cmath>
#include <memory>

using std::shared_ptr;
using std::vector;

/// Generic kernel class
template< typename _T>
class CKernel
{
public:
  /// Constructor
  explicit CKernel( bool normalized = false ) {bNorm = normalized;}

  /// Destructor
  virtual ~CKernel() {}

private:
  /// A flag indicating whether to normalize the kernel output from CKernel::eval
  bool bNorm;

public:
  /// Takes two objects and returns the kernel value applying any necessary post-processing (such as normalization)
  virtual double operator() ( const _T& pobj1, 
			      const _T& pobj2 ) const
  {
    if( bNorm == false ) return eval( pobj1, pobj2 );
    else
      {
	double n = eval( pobj1, pobj2 );
	double d = eval( pobj1, pobj1 ) * eval( pobj2, pobj2 );
	if( d == 0.0 ) return n;
	else return n / std::sqrt(d);
      }
  }

private:
  /// Takes two objects and returns the kernel value with no post-processing. Must be overridden by children classes
  virtual double eval( const _T& pobj1,
		       const _T& pobj2 ) const = 0;

};

////////////////////////////////////////////////////////////////////////////////////

/// A simple identity kernel that returns 1 if the objects are equal and 0 otherwise
template< typename _T >
class CIdentityKernel : public CKernel<_T>
{
public:
  /// Constructor
  CIdentityKernel() : CKernel<_T>(false) {}

private:
  /// Takes two objects and returns the kernel value with no post-processing. Must be overridden by children classes
  virtual double eval( const _T& pobj1,
		       const _T& pobj2 ) const
  {
    if( pobj1 == pobj2 ) return 1.0;
    else return 0.0;
  }
  
};

////////////////////////////////////////////////////////////////////////////////////

/// Gaussian kernel class that can be used on top of another kernel
template< typename _T >
class CGaussianKernel : public CKernel< _T >
{
public:
  /// Constructor
  CGaussianKernel( shared_ptr< CKernel<_T> const > k,
		   double gamma_param = 0.1,
		   bool normalized = false ) : CKernel< _T >(normalized),
					       K(k)
  {
    gamma = gamma_param;
  }

private:
  /// The width parameter
  double gamma;

  /// The underlying kernel used to carry out internal computations
  shared_ptr< CKernel<_T> const > K;
  
private:
  /// Takes pointers to two objects and returns the kernel value
  virtual double eval( const _T& pobj1,
		       const _T& pobj2 ) const
  {
    double dist_sq = (*K)(pobj1, pobj1) - 2.0 * (*K)(pobj1, pobj2) + (*K)(pobj2, pobj2);
    return std::exp( -gamma * dist_sq );
  } 
};

////////////////////////////////////////////////////////////////////////////////////

/// A generic composite kernel defined over multiple feature spaces
template< typename _T >
class CCompositeKernel : public CKernel< vector< shared_ptr<_T> > >
{
private:
  typedef vector< shared_ptr<_T> > _U;
  typedef shared_ptr< const CKernel< _T > > _pK;

public:
  /// Default constructor
  CCompositeKernel()
    : CKernel<_U>(false) {}
  
  /// Constructor
  explicit CCompositeKernel( const vector< _pK >& k  )
    : CKernel<_U>(false),
      kernels( k ) {}

  /// Destructor
  virtual ~CCompositeKernel() {}

private:
  /// A collection of kernels associated with the feature spaces
  vector< _pK > kernels;

public:
  /// Clears the associated vector of kernels
  void clear() {kernels.clear();}

  /// Appends another feature space kernel
  void addKernel( const _pK& k ) { kernels.push_back(k); }

  /// Returns the number of internal kernels
  unsigned int size() const {return kernels.size();}

private:
  /// Takes two samples that span all feature spaces and computes the composite kernel value
  virtual double eval( const _U& pobj1,
		       const _U& pobj2 ) const
  {
    if( pobj1.size() != pobj2.size() || pobj1.size() != kernels.size() )
      throw std::logic_error( "Mismatched dimensionality in CCompositeKernel" );
    
    double res = 0.0;
    for( unsigned int i = 0; i < kernels.size(); i++ )
      res += (*kernels[i])( *(pobj1[i]), *(pobj2[i]) );
    return res;
  }
};

////////////////////////////////////////////////////////////////////////////////////

/// A product joint kernel of the form KJ = KX * KY
class CProdJointKernel
{
public:
  /// Constructor
  explicit CProdJointKernel() {}

public:
  /// Takes a value of the input kernel KX and the output kernel KY and returns the joint kernel value
  virtual double operator() ( double KX, double KY ) const {return KX * KY;}
};

////////////////////////////////////////////////////////////////////////////////////

/// A polynomial joint kernel of the form KJ = (KX + KY + 1)^d
class CPolyJointKernel
{
public:
  /// Constructor
  explicit CPolyJointKernel( unsigned int degree ) : d(degree) {}

private:
  /// Degree of the polynomial
  unsigned int d;

public:
  /// Takes a value of the input kernel KX and the output kernel KY and returns the joint kernel value
  virtual double operator() ( double KX, double KY ) const {return std::pow( (KX + KY + 1.0), static_cast<double>(d) );}
};

////////////////////////////////////////////////////////////////////////////////////

/// A homogeneous polynomial joint kernel of the form KJ = (KX + KY)^d
class CPolyHomJointKernel
{
public:
  /// Constructor
  explicit CPolyHomJointKernel( unsigned int degree ) : d(degree) {}

private:
  /// Degree of the polynomial
  unsigned int d;

public:
  /// Takes a value of the input kernel KX and the output kernel KY and returns the joint kernel value
  virtual double operator() ( double KX, double KY ) const {return std::pow( (KX + KY), static_cast<double>(d) );}
};

#endif
