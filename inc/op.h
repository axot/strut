// -*-c++-*-
/// \file op.h
/// \brief interface to the operator classes
/// \author Artem Sokolov

#ifndef OP_H__INCLUDED
#define OP_H__INCLUDED

#include <memory>
using std::shared_ptr;

/// Generic less-than operator
template< typename _T>
class CLTOp
{
public:
  /// Constructor
  CLTOp() {}

  /// Destructor
  virtual ~CLTOp() {}

public:
  /// Returns true if a is "strictly less" than b
  virtual bool operator() ( shared_ptr< const _T > a,
			    shared_ptr< const _T > b ) const = 0;
};

/// Generic equal operator
template< typename _T >
class CEQOp
{
public:
  /// Constructor
  CEQOp() {}

  /// Destructor
  virtual ~CEQOp() {}

public:
  /// Returns true if a "equals" b
  virtual bool operator() ( shared_ptr< const _T > a,
			    shared_ptr< const _T > b ) const = 0;
};

/// Uses the inherent operator< of the objects
template< typename _T >
class CInherentLTOp : public CLTOp< _T >
{
public:
  /// Constructor
  CInherentLTOp() {}

public:
  /// Returns true if a is "strictly less" than b
  virtual bool operator() ( shared_ptr< const _T > a,
			    shared_ptr< const _T > b ) const
  {
    return ((*a) < (*b));
  }

};

/// Uses the inherent operator== of the objects
template< typename _T >
class CInherentEQOp : public CEQOp< _T >
{
public:
  /// Constructor
  CInherentEQOp() {}

public:
  /// Returns true if a is "equal" to b
  virtual bool operator() ( shared_ptr< const _T > a,
			    shared_ptr< const _T > b ) const
  {
    return ((*a) == (*b));
  }
};

#endif
