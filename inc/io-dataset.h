// -*-c++-*-
/// \file io-dataset.h
/// \brief Interface to CIODataSet
/// \author Artem Sokolov

#ifndef IO_DATASET_H__INCLUDED
#define IO_DATASET_H__INCLUDED

#include "dataset.h"

#include <memory>
#include <functional>
using std::shared_ptr;
using std::function;

/// An implementation of pairing of objects from arbitrary input and output spaces
template< typename _I, typename _O >
class CIODataSet
{
public:
  /// Constructor
  CIODataSet( const function< double( const _I&, const _I& ) >& fkeri,
	      const function< double( const _O&, const _O& ) >& fkero,
	      const function< double( const _O&, const _O& ) >& floss,
	      const function< double( double, double ) >& f_ioker );
  
  /// Constructor
  CIODataSet( shared_ptr< CDataSet<_I> > pi,
  	      shared_ptr< CDataSet<_O> > po,
	      const vector< unsigned int >& io_map,
	      const function< double( double, double ) >& f_ioker );

  /// Constructor
  CIODataSet( const function< double( const _I&, const _I& ) >& fkeri, 
	      shared_ptr< CDataSet<_O> > po,
	      const function< double( double, double ) >& f_ioker );

  /// Copy constructor
  CIODataSet( const CIODataSet& other );

private:
  /// Default constructor
  CIODataSet( const function< double( double, double ) >& f_ioker )
    : fioker( f_ioker ) {}

public:
  /// Destructor
  ~CIODataSet() {}

private:
  /// The associated dataset over the input space
  shared_ptr< CDataSet<_I> > pids;

  /// The associated dataset over the output space
  shared_ptr< CDataSet<_O> > pods;

  /// Maps samples of pids -> samples of pods (many-to-one mapping)
  vector< unsigned int > iomap;

  /// The associated joint kernel
  const function< double( double, double ) > fioker;

public:
  /// Direct, but const access to CIODataSet::pids
  const shared_ptr< const CDataSet<_I> > getI() const {return pids;}

  /// Direct, but const access to CIODataSet::pods
  const shared_ptr< const CDataSet<_O> > getO() const {return pods;}

public:
  /// Splits the dataset into training and test subsets
  void splitTrainTest( virange_t vTrain, virange_t vTest,
		       shared_ptr< CIODataSet<_I, _O> >& pTrain,
		       shared_ptr< CIODataSet<_I, _O> >& pTest ) const;

  /// Returns the number of data samples in the set
  unsigned int sizeI() const {return pids->size();}

  /// Returns the number of unique examples in the output space
  unsigned int sizeO() const {return pods->size();}

  /// Subsamples the specified samples of the dataset
  void subsample( const vector< unsigned int >& indices );

  /// Randomly shuffles the examples
  void random_shuffle();

  /// Matches up the sample IDs between i_ds and o_ds and stores the matches internally
  /** \param[in] i_ds a dataset over the input space
      \param[in] o_ds a dataset over the output space
   */
  unsigned int addSets( shared_ptr<CDataSet<_I> const> i_ds, 
			shared_ptr<CDataSet<_O> const> o_ds )
  { return addSets( *i_ds, *o_ds ); }

  /// Matches up the sample IDs between i_ds and o_ds and stores the matches internally
  /** \param[in] i_ds a dataset over the input space
      \param[in] o_ds a dataset over the output space
   */
  unsigned int addSets( const CDataSet<_I>& i_ds, 
			const CDataSet<_O>& o_ds );

  /// Add an input-output sample pairing
  void addSample( const string& name, shared_ptr<_I const> pISample,
		  shared_ptr<_O const> pOSample );

  /// Adds a sample to the input space
  void addInputSample( const string& name, shared_ptr< _I const > pSample,
		       unsigned int mapping );

  /// Adds a sample to the output space, maintaining uniqueness
  int addOutputSample( const string& name, shared_ptr<_O const> pOSample );

  /// Maps an index into the input dataset to an index into the output dataset
  unsigned int map( unsigned int i ) const {return iomap[i];}

  /// Remaps the i^th input sample to k^th output sample
  void remap( unsigned int i, unsigned int k ) {iomap[i] = k;}

  /// Recaches the kernel matrices
  void cache();

  /// Caches the input space against an external dataset
  void cacheIExternal( shared_ptr<CDataSet<_I> const > piExt ) const
  { pids->cacheExternal( piExt ); }

  /// Looks up the sample index by the sample ID, returns -1 if no such sample
  int s2i( const string& s ) const {return pids->s2i(s);}

  /// Looks up the sample ID by its index
  const string i2s( unsigned int i ) const {return pids->i2s(i);}

  /// Computes the input-space kernel value
  double ikernel( unsigned int xi, unsigned int xj ) const
  { return pids->kernel( xi, *pids, xj ); }

  /// Computes the input-space kernel value
  double ikernel( unsigned int xi, std::pair< const CDataSet<_I>&, unsigned int > x ) const
  { return pids->kernel( xi, x.first, x.second ); }

  /// Computes the input-space kernel value
  double ikernel( unsigned int xi, const CDataSet<_I>& other, unsigned int xj ) const
  { return pids->kernel( xi, other, xj ); }

  /// Computes the output-space kernel value
  double okernel( unsigned int yi, unsigned int yj ) const
  { return pods->kernel( yi, *pods, yj ); }

  /// Computes the joint kernel using pre-computed kernel values for the input and output spaces
  double iokernel( double xker, double yker ) const
  { return fioker( xker, yker ); }

  /// Computes the joint kernel using the associated object
  template< typename _TX, typename _TY >
  double iokernel( unsigned int xi1, unsigned int yi1, _TX xi2, _TY yi2 ) const
  { return fioker( ikernel( xi1, xi2 ), okernel( yi1, yi2 ) ); }

  /// Computes the loss value between output-space samples i and j, or retrieves a cached value
  double oloss( unsigned int yi, unsigned int yj ) const
  { return pods->loss( yi, yj ); }

  /// Computes the loss value between an output-space sample i and a foreign sample
  double oloss( unsigned int yi, shared_ptr< _O const > py ) const
  { return pods->loss( yi, py ); }

};

#include "io-dataset_impl"

#endif

