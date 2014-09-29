// -*-c++-*-
/// \file infer.h
/// \brief Inference algorithms
/// \author Artem Sokolov

#ifndef INFER_H__INCLUDED
#define INFER_H__INCLUDED

#include "types.h"

// Forward declarations
template< typename _I, typename _O > class CClassifier;
template< typename _I, typename _O > class CBLASTNN;
template< typename _I, typename _O > class C1NN;
template< typename _I, typename _O, char _R > class CnsSSVM;

/// Generic inference that loops over all available labels
template< typename _I, typename _O >
struct FLoopyArgmax
{
  /// The associated classifier
  const CClassifier<_I, _O>& clsf;

  /// Constructor
  FLoopyArgmax( const CClassifier<_I, _O>& clsf_ ) : clsf( clsf_ ) {}

  /// Infers the most compatible label for xi'th sample from dataset ds
  unsigned int operator()( pair< const CDataSet<_I>&, unsigned int > smpl );
};

/// Infernce for the BLAST-NN classifier
template< typename _I, typename _O >
struct FBLASTArgmax
{
  /// The associated classifier
  const CBLASTNN<_I,_O>& clsf;

  /// Constructor
  FBLASTArgmax( const CBLASTNN<_I,_O>& clsf_ ) : clsf( clsf_ ) {}

  /// Infers the most compatible label for xi'th sample from dataset ds
  unsigned int operator()( pair< const CDataSet<_I>&, unsigned int > smpl );
};

/// Inference over a GO hierarchy using a linear kernel model
template< typename _I, char _R >
struct SGOArgmax
{
  /// The associated classifier
  const CnsSSVM<_I,CSparseSample,_R>& svm;
  
  /// Constructor
  SGOArgmax( const CnsSSVM<_I,CSparseSample,_R>& svm_ ) : svm( svm_ ) {}

  /// Infers the most compatible label for xi'th sample from dataset ds
  unsigned int operator()( pair< const CDataSet<_I>&, unsigned int > smpl );
};

#include "infer_impl"

#endif
