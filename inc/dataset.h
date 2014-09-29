// -*-c++-*-
/// \file dataset.h
/// \brief Interface to the generic CDataSet
/// \author Artem Sokolov

#ifndef DATASET_H__INCLUDED
#define DATASET_H__INCLUDED

#include "kernel.h"
#include "loss.h"
#include "op.h"
#include "types.h"
#include "misc.h"

#include <memory>
#include <string>

using std::auto_ptr;
using std::shared_ptr;
using std::string;
using std::vector;

/// A generic dataset container, completely decoupled from the feature space
template< typename _T >
class CDataSet
{
public:
  typedef function< double( const _T&, const _T& ) > binop_t;

private:
  /// Invalidates cache of the dataset when the object goes out of scope
  class CCacheInvalidator
  {
  public:
    /// Constructor
    CCacheInvalidator( CDataSet<_T>* ds ) : pds( ds ) {}

    /// Destructor
    ~CCacheInvalidator() { pds->invalidateCache(); }

  private:
    /// The associated dataset
    CDataSet<_T>* const pds;
  };

public:
  /// Constructor
  CDataSet( const binop_t& f_ker = CIdentityKernel<_T>(),
	    const binop_t& f_loss = CIdentityLoss<_T>() )
    : fker( f_ker ), floss( f_loss ), isCached( false ) {}

  /// Copy constructor
  CDataSet( const CDataSet& other );

  /// Assignment operator
  CDataSet& operator=( const CDataSet& other );

  /// Non-virtual destructor to discourage deriving classes
  ~CDataSet() {}

private:
  /// Maps indices to sample IDs
  vector< string > sampleIDs;		
  
  /// Maps sample IDs to indices
  simap_t i_sampleIDs;		

  /// The samples
  vector< shared_ptr<_T> > samples;

  /// A pointer to the associated CKernel object
  binop_t fker;

  /// A pointer to the associated CLoss object
  binop_t floss;

  /// Signifies whether the current cache is valid
  bool isCached;

  /// A cache of the kernel matrix
  vector< double > kernelCache;

  /// A cache of the loss matrix
  vector< double > lossCache;

  /// A pointer to the dataset associated with the external cache
  shared_ptr< CDataSet<_T> const > pExt;

  /// A cache of kernel values against an external dataset, indexed as (local index, external index)
  vector< double > externalCache;

public:
  /// Displays the dataset to an arbitrary output stream
  void display( std::ostream& os = std::cout ) const
  { for( unsigned int i = 0; i < samples.size(); i++ ) displaySample( i, os ); }

  /// Displays local samples to an arbitrary output stream
  void displaySamples( vector< string > ids, std::ostream& os = std::cout ) const;

  /// Displays the i^th sample to an arbitrary output stream
  void displaySample( unsigned int i, std::ostream& os = std::cout ) const
  { os<<i2s(i)<<","<<*samples[i]; }

  /// Saves the dataset to a file
  void save( const string& filename ) const
  { std::ofstream ofs( filename.c_str() ); display(ofs); }

public:
  /// Returns a shared pointer to a const sample at the i^th position
  const shared_ptr<_T const> getSample( unsigned int i ) const {return samples[i];}

  /// Returns a shared pointer to a non-const sample at the i^th position
  const shared_ptr<_T> getSampleMod( unsigned int i ) {invalidateCache(); return samples[i];}

  /// Returns a shared pointer to the associated kernel
  const binop_t getKernel() const {return fker;}

  /// Returns the associated loss
  const binop_t getLoss() const {return floss;}

public:
  /// Clears the dataset
  virtual void clear();

  /// Adds a new sample to the data set
  /**
     The ownership of the sample is transfered to the dataset
     Use CDataSet::addSampleConst to make a copy of the sample if the original must be preserved
     \param[in] sampleID the ID of the new sample
     \param[in] smpl a pointer to the actual sample
     \param[in] bOverwrite specifies whether a sample with sampleID is to be overwritten if it already exists
  */
  void addSample( const string& sampleID,
		  auto_ptr<_T> smpl,
		  bool bOverwrite = false )
  {
    shared_ptr<_T> p( smpl );
    addSample( sampleID, p, bOverwrite );
  }

private:
  /// Adds a new pre-created sample to the data set
  void addSample( const string& sampleID,
  		  shared_ptr<_T> smpl,
		  bool bOverwrite = false );

public:
  /// Adds a new sample to the data set
  /**
     The function first makes a copy of the provided sample, then adds the copy to the dataset.
     The original sample is preserved
     \param[in] sampleID the ID of the new sample
     \param[in] smpl A pointer to the actual sample. A copy of this sample will be added to the dataset
     \param[in] bOverwrite specifies whether a sample with sampleID is to be overwritten if it already exists
  */
  void addSampleConst( const string& sampleID,
		       shared_ptr<const _T> smpl,
		       bool bOverwrite = false );

  /// Adds the i^th sample of the other dataset to this one
  void addSample( shared_ptr< CDataSet<_T> const > other,
  		  unsigned int i, bool bOverwrite = false,
  		  const string& sampleID = "" );

  /// Adds another dataset to this one
  void addSet( shared_ptr< CDataSet<_T> const > other,
	       bool bOverwrite = false,
	       const string& appendID = "" );

public:
  /// Returns the size of the dataset
  unsigned int size() const { return samples.size(); }

  /// Subsamples the dataset by indices
  void subsample( const vector< unsigned int >& indices );

  /// Subsamples the dataset by sample IDs
  /**
     \param[in] strict if false, missing ids are ignored; exception is thrown, if true
   */
  void subsample( const vector< string >& ids, bool strict = true );

public:
  /// Returns all sample IDs
  const vector< string > getSampleIDs() const { return sampleIDs; }

  /// Looks up the sample index by the sample ID, returns -1 if no such sample
  int s2i( const string& ) const;

  /// Looks up the sample ID by its index
  const string i2s( unsigned int i ) const;

  /// Renames the i^th sample to s
  void rename( unsigned int i, const string& s );

  /// Returns index of the sample or -1 if there is no match
  int findSample( shared_ptr< _T const > smpl ) const;

public:
  /// Computes the kernel value between the i^th sample and a sample from another (or the same) dataset
  double kernel( unsigned int i, const CDataSet<_T>& other, unsigned int j ) const;

  /// Computes the loss between the i^th sample and the second argument
  double loss( unsigned int i, shared_ptr< _T const > smpl ) const
  { return floss( *(samples[i]), *smpl ); }

  /// Computes the loss between the i^th and j^th samples
  double loss( unsigned int i, unsigned int j ) const;

  /// Computes the average loss per sample, given the truth
  double loss( const CDataSet<_T>& truth ) const;

  /// Caches the kernel and the loss matrices
  void cache();

  /// Caches the kernel values with an external dataset
  void cacheExternal( shared_ptr< const CDataSet<_T> > pExtDS );

private:
  /// Clears the current cache
  void invalidateCache()
  { isCached = false; 
    pExt = shared_ptr< CDataSet<_T> >( static_cast< CDataSet<_T>* >(NULL) ); }
};

/// Computes a set of common sample IDs with another dataset
template< typename _T >
vector< string > commonSampleIDs( const CDataSet<_T>& pds1,  
				  const CDataSet<_T>& pds2 );

#include "dataset_impl"

#endif
