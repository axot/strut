// -*-c++-*-
/// \file loss.h
/// \brief Interface to loss function objects
/// \author Artem Sokolov

#ifndef LOSS_H__INCLUDED
#define LOSS_H__INCLUDED

#include "kernel.h"

#include <memory>
using std::shared_ptr;

/// Generic loss function object
template< typename _T >
class CLoss
{
public:
  /// Constructor
  CLoss() {}

  /// Destructor
  virtual ~CLoss() {}

public:
  /// Takes two objects and returns the loss value
  virtual double operator() ( const _T& pobj1, 
			      const _T& pobj2 ) const = 0;
};

/// Generic kernel loss
template< typename _T >
class CKernelLoss : public CLoss<_T>
{
public:
  typedef function< double( const _T&, const _T& ) > binop_t;

public:
  /// Constructor
  CKernelLoss( const binop_t& fk ) : fker( fk ) {}

private:
  /// The associated kernel
  binop_t fker;

public:
  /// Takes two objects and returns the loss value
  virtual double operator() ( const _T& pobj1, 
			      const _T& pobj2 ) const
  {
    double kn = fker(pobj1, pobj2);
    double k1 = fker(pobj1, pobj1);
    double k2 = fker(pobj2, pobj2);
    if( (k1 + k2) == 0.0 ) return 1.0;
    else return (1.0 - 2.0 * kn / (k1 + k2));
  } 
  
};

/// A simple identity loss that returns 0 if the objects are equal and 1 otherwise
template< typename _T >
class CIdentityLoss : public CLoss<_T>
{
public:
  /// Constructor
  CIdentityLoss() {}

public:
  /// Takes two objects and returns the loss value
  virtual double operator() ( const _T& obj1, 
			      const _T& obj2 ) const
  {
    if( obj1 == obj2 ) return 0.0;
    else return 1.0;
  }
};


#endif
