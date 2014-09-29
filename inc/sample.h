// -*-c++-*-
/// \file sample.h
/// \brief headers for the Sample Classes
/// \author Artem Sokolov

#ifndef SAMPLE_H__INCLUDED
#define SAMPLE_H__INCLUDED

#include "types.h"
#include "featmap.h"
#include "dataset.h"
#include "kernel.h"

/// Sparse sample representation
/** The internal representation is a collection of (key, value) pairs
    Class invariant: the sample must be sorted by the keys
    Class invariant: zero values are not represented
    Class invariant: a single value is associated with each key
*/
class CSparseSample
{
private:

  // The associated feature map
  shared_ptr< CFeatMap const > pfmap;
  
  /// (key, value) pairs
  std::vector< std::pair< unsigned int, double > > content;

public:
  /// Constructor
  explicit CSparseSample( shared_ptr< CFeatMap const > fmap );

  /// Displays the sample to an arbitrary output stream
  void display( std::ostream& os = std::cout ) const;

  /// Returns the associated feature map
  const shared_ptr< CFeatMap const > getFeatMap() const {return pfmap;}

public:
  /// Returns the value for the requested key
  double getValue( string key ) const;

  /// Returns the value for the requested key
  double getValue( unsigned int key ) const;

  /// Sets the value for the requested key
  void setValue( unsigned int key, double val );

  /// Sets the value for the requested key
  void setValue( string key, double val );
  
  /// Maps the sample to a new feature map
  void remap( shared_ptr< CFeatMap const > newmap );

public:
  /// Computes the dot product with another sparse sample
  double dot( const CSparseSample& other ) const;

  /// Homogeneous version of the dot product
  /** This version subtracts 1 from the dot product, which is useful
      if the two samples come from the same ontology, for example,
      and where the fact that they share a root node is irrelevant
   */
  double dot_hom( const CSparseSample& other ) const
  { return dot( other ) - 1.0; }

  /// Multiplies the sparse vector by a scalar
  CSparseSample& operator*=( const double& coeff );

  /// Adds another vector to *this
  CSparseSample& operator+=( const CSparseSample& rhs );

  /// Multiplies the sparse vector by a scalar
  const CSparseSample operator*( const double& coeff ) const;

  /// Lexicographic key-based less-than operator
  bool operator<( const CSparseSample& other ) const;

  /// Lexicographic key-based equality operator
  bool operator==( const CSparseSample& other ) const;

  /// Returns the number of non-zero entries in the sample
  unsigned int L0() const { return content.size(); }
};

/// A linear kernel over CSparseSample objects
class CSparseKernel : public CKernel< CSparseSample >
{
public:
  /// Constructor
  explicit CSparseKernel( bool normalized ) : CKernel< CSparseSample >(normalized) {}

private:
  /// Takes pointers to two objects and returns the kernel value
  virtual double eval( const CSparseSample& pobj1,
		       const CSparseSample& pobj2 ) const
  { return pobj1.dot( pobj2 ); }

};

/// A linear homogeneous kernel over CSparseSample objects
class CSparseHomKernel : public CKernel< CSparseSample >
{
public:
  /// Constructor
  explicit CSparseHomKernel( bool normalized ) : CKernel< CSparseSample >(normalized) {}

private:
  /// Takes two objects and returns the kernel value
  virtual double eval( const CSparseSample& pobj1,
		       const CSparseSample& pobj2 ) const
  { return pobj1.dot_hom( pobj2 ); }

};

/// I/O
std::ostream& operator<< ( std::ostream& os, const CSparseSample& smpl );

////////////////////////// Values //////////////////////////////

/// Returns the range of values in the dataset
std::pair< double, double > getRange( const CDataSet<CSparseSample>& ds );

////////////////////////// Features ////////////////////////////

/// Returns the number of features in a dataset of sparse samples
unsigned int nFeats( const CDataSet<CSparseSample>& ds );

/// Returns the number of features in a dataset of sparse samples
unsigned int nFeats( const shared_ptr< const CDataSet<CSparseSample> >& pds );

/// Returns the number of samples that have feature f
unsigned int nSamplesWFeat( const CDataSet<CSparseSample>& ds, string f );

/// Returns a mapping of features to the number of samples that have it
void nSamplesWFeat( const CDataSet<CSparseSample>& ds, shared_ptr< CFeatMap const > pfmap, simap_t& counts );

/// Computes a set of common feature IDs between two datasets
vector< string > commonFeatIDs( const CDataSet<CSparseSample>& ds1,
				const CDataSet<CSparseSample>& ds2 );

/// Remaps the dataset to a new feature map
void remap( CDataSet<CSparseSample>& ds, shared_ptr< CFeatMap > fmap );

/////////////////////////// Manipulation ////////////////////////

/// Removes all samples that have fewer non-zero features than thresh
void cropSamples( unsigned int thresh, CDataSet<CSparseSample>& ds );

/// Thresholds the dataset, setting the top-scoreing k samples for each feature to 1 and rest to 0
/**
   \param[in] A thresholding profile that contains (feature ID, number of top-scoring samples) pairs
*/
void threshold( CDataSet<CSparseSample>& ds, std::map< string, unsigned int > prof );

/// Thresholds the dataset, setting all entries for feature k that (strictly) exceed prof[k] to 1 and rest to 0
void threshold( CDataSet<CSparseSample>& ds, std::map< string, double > prof );

/// Adds new features to the samples of ds, overwriting any previous values
void expand( CDataSet<CSparseSample>& ds, shared_ptr< CFeatMap > pfm,
	     const CDataSet<CSparseSample>& other );

////////////////////////// Statistics ////////////////////////////

#include "clsf.h"

/// Computes the average precision and recall per feature, given the truth
std::pair< double, double > computePnR( const CDataSet<CSparseSample>& dsPred, const CDataSet<CSparseSample>& dsTruth );

/// Computes prediction scores for each output-space features. Save the results to fnOutput
template< typename _I >
void predScores( const CClassifier<_I, CSparseSample>& clsf, shared_ptr< const CDataSet<_I> > pdsTest, const string& fnOutput );

/// Computes prediction scores for each output-space features. Returns the results as a sparse dataset
template< typename _I >
CDataSet<CSparseSample> predScores( const CClassifier<_I, CSparseSample>& clsf, shared_ptr< const CDataSet<_I> > pdsTest );

////////////////// Multiple Kernel Learning ////////////////////

typedef vector< CSparseSample > vSparseSample;

/// A composite kernel that is the sum of individual kernels, each of which is normalized
class CCompositeSparseKernel : public CKernel< vSparseSample >
{
public:
  /// Constructor
  explicit CCompositeSparseKernel( bool normalized ) : CKernel< vSparseSample >(normalized), ker( true ) {}

private:
  /// Kernel used for individual spaces
  CSparseKernel ker;

private:
  /// Takes two objects and returns the kernel value
  virtual double eval( const vSparseSample& obj1,
		       const vSparseSample& obj2 ) const
  {
    if( obj1.size() != obj2.size() )
      throw std::logic_error( "Incompatible computation in CCompositeSparseKernel: different number of kernels available for each sample" );

    // Traverse the kernel spaces
    double res = 0.0;
    for( unsigned int i = 0; i < obj1.size(); ++i )
      res += ker( obj1[i], obj2[i] );
    return res;
  }
  
};

/// Adds new features to the samples of ds, appending them as another feature space
void expand( CDataSet<vSparseSample>& ds, const CDataSet<CSparseSample>& other, bool bRemoveMissing = false );

/// Returns the number of features in an MKL dataset
unsigned int nFeats( const CDataSet<vSparseSample>& ds );

/// Returns the number of kernels used by an MKL dataset
unsigned int nKernels( const CDataSet<vSparseSample>& ds );

///////////////////////// Gene Ontology //////////////////////////

#include "go-annotation.h"
#include "blastout.h"

/// Takes a set of protein names and creates a sparse dataset over the annotations for those proteins from a CGOACollection
void makeSparseDataset( const GO::CGOACollection& goaSource,
		        const vector<string>& protnames,
			CDataSet<CSparseSample>& ds,
			const GO::CGOContainer& goGraph,
			GO::ONTOLOGY_INDEX filter = (GO::GO_MF | GO::GO_BP | GO::GO_CC),
			shared_ptr< CFeatMap > pfmap = 
			shared_ptr< CFeatMap >( static_cast<CFeatMap*>(NULL)) );

/// Takes a set of protein names and creates a sparse dataset over the BLAST hits for those proteins
/** Uses negative log e-values as features
    \param[in] source An object containing the BLAST hits
    \param[in] protnames A set of proteins over which the dataset is created
    \param[out] ds The result is appended to this dataset  
    \param[in] lower_thresh All e-values below this value are considered to be equal to it
    \param[in] upper_thresh All e-values above this value are ignored
    \param[in] bAddEmpty If set to true, then empty entries are added for proteins that are not found in source
    \param[in] pfmap A pointer to the feature map to use. New map is created if NULL.
    */
void makeSparseDataset( const GO::CBLASTOutput& source,
			CDataSet<CSparseSample>& ds,
			double lower_thresh, double upper_thresh,
			shared_ptr< CFeatMap > pfmap = 
			shared_ptr< CFeatMap >( static_cast<CFeatMap*>(NULL)) );

#include "sample_impl"

#endif
