// -*-c++-*-
/// \file blast-nn.h
/// \brief BLAST-nearest neighbor
/// \author Artem Sokolov

#ifndef BLAST_NN_H__INCLUDED
#define BLAST_NN_H__INCLUDED

#include "blastout.h"
#include "clsf.h"

#include <memory>

using std::shared_ptr;

/// BLAST-based nearest neighbor
template< typename _I, typename _O >
class CBLASTNN : public CClassifier<_I, _O >
{
public:
  /// Constructor
  explicit CBLASTNN( shared_ptr< GO::CBLASTOutput const > p_bo )
    : CClassifier<_I, _O>( "BLAST-NN", FBLASTArgmax<_I,_O>(*this) ),
      pbo( p_bo ) {}

private:
  /// The associated CBLASTOutput object
  shared_ptr< GO::CBLASTOutput const > pbo;

private:
  /// Trains the classifier  
  virtual void train() {}

public:
  /// Clear the learned information for the classifier
  virtual void clear() {}  

  /// Generic implementation of the compatibility function is not used
  virtual double f( pair< const CDataSet<_I>&, unsigned int > x, unsigned int dso ) const
  { throw std::logic_error( "Function is not used" ); }

public:
  /// Finds the closest kernel-based neighbor of sample dsi
  /** \return index into the input-space dataset of the closest neightbor, or -1 if no match
   */
  template< typename _T >
  int findKernelNeighbor( _T x ) const;

  /// Finds the closest BLAST neightbor in the training data for a particular protein
  /** \return index into the input-space dataset of the closest neightbor, or -1 if no match
   */
  int findBLASTNeighbor( string strName ) const;

  /// Returns the name of an external sample
  std::string findName( std::pair< const CDataSet<_I>&, unsigned int > dsx ) const
  { return dsx.first.i2s( dsx.second ); }
};

template< typename _I, typename _O >
int CBLASTNN<_I,_O>::findBLASTNeighbor( string strName ) const
{
  int res = -1;

  // Find the associated BLAST hits
  GO::CBLASTOutput::emap_const_iter_t iter = pbo->find( strName );
  if( iter != pbo->end() )
    {
      // Find the closest hit
      double best_e_val = std::numeric_limits<double>::infinity();
      for( unsigned int j = 0; j < iter->second.size(); j++ )
	{
	  // Retrieve the entry
	  string jStr = iter->second[j].subject_id;
	  double jVal = iter->second[j].e_value;

	  // Ensure presence of the name in the dataset
	  int dsj = this->pdsTrain->s2i( jStr );
	  if( dsj < 0 ) continue;

	  // Compare to the best found
	  if( jVal < best_e_val )
	    {
	      best_e_val = jVal;
	      res = dsj;
	    }
	}
    }

  return res;
}

template< typename _I, typename _O >
template< typename _T >
int CBLASTNN<_I,_O>::findKernelNeighbor( _T x ) const
{
  int res = -1;

  double best_K = -1.0 * std::numeric_limits<double>::infinity();

  // Infer the closest sample in the input space (1-NN)
  for( unsigned int i = 0; i < this->pdsTrain->sizeI(); i++ )
    {
      // Compute the kernel between the training and test samples
      double K = this->pdsTrain->ikernel( i, x );

      // Compute argmax
      if( K > best_K )
	{
	  best_K = K;
	  res = i;
	}
    }

  return res;
}

#endif
