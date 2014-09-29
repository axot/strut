// -*-c++-*-
/// \file blast-wrap.h
/// \brief interface to BLAST-based nearest-neighbor
/// \author Artem Sokolov

#ifndef BLAST_WRAP_H__INCLUDED
#define BLAST_WRAP_H__INCLUDED

#include "blastout.h"
#include "io-dataset.h"

/// Generates a set of output sample indices for a given input space index
class COutputSubspaceGenerator
{
protected:
  typedef const vector<unsigned int>::const_iterator iter_t;

public:
  /// Constructor
  explicit COutputSubspaceGenerator() {}
  
  /// Destructor
  virtual ~COutputSubspaceGenerator() {}

private:
  /// The set of output labels which manifest the subspace. May be changed by operator()
  vector<unsigned int> indices;

protected:
  /// Adds a new index to the subspace
  void add( unsigned int i ) {indices.push_back(i);}

  /// Returns an iterator to the beginning
  iter_t begin() const {return indices.begin();}

  /// Returns an iterator to the end
  iter_t end() const {return indices.end();}

private:
  /// Generates the output subspace for a given input space sample id
  virtual void generate( string sample_id ) = 0;

public:
  /// Returns a (begin, end) pair of output-space iterators for a given input sample id
  virtual std::pair< iter_t, iter_t > operator() (string sample_id)
  { 
    indices.clear();
    generate( sample_id );
    return std::pair<iter_t, iter_t>( indices.begin(), indices.end() ); 
  }
};


template< typename _I, typename _O >
class CBLASTOSG : public COutputSubspaceGenerator
{
private:
  typedef const vector< unsigned int >::const_iterator iter_t;

public:
  /// Constructor
  explicit CBLASTOSG( shared_ptr< const CIODataSet<_I, _O> > p_ds,
		      shared_ptr< GO::CBLASTOutput const > p_bo,
		      GO::CGOContainer& goCont,
		      bool bRecombine = false )
    : pds( p_ds ), pbo( p_bo ), goGraph( goCont ), recomb( bRecombine )
  {}

  /// Destructor
  virtual ~CBLASTOSG() {}

private:
  /// The associated dataset
  shared_ptr< const CIODataSet<_I, _O> > pds;

  /// The associated CBLASTOutput object
  shared_ptr< const GO::CBLASTOutput > pbo;

  /// The associated GO graph
  const GO::CGOContainer& goGraph;

  /// Whether of not the outputs should be recombined
  bool recomb;

private:
  /// Generates the output subspace for a given input space index
  virtual void generate( string sample_id );

  /// Recombines the BLAST outputs
  void recombineOutputs( const std::set< unsigned int >& po );

};

#include "blast-wrap_impl"


#endif
