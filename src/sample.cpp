// sample.cpp - implementation of Sample Classes
//
// by Artem Sokolov

#include "sample.h"
#include <limits>

// Construction
CSparseSample::CSparseSample( shared_ptr< CFeatMap const > fmap )
  : pfmap( fmap ) { if( !pfmap ) throw std::logic_error( "No feature map provided to a CSparseSample constructor" ); }

// Displays the sample to an arbitrary output stream
void CSparseSample::display( std::ostream& os ) const
{
  for( unsigned int i = 0; i < content.size(); i++ )
    {
      string s = pfmap->i2f( content[i].first );
      if( s.empty() ) throw std::logic_error( "Unmappable feature in CSparseSample::display" );
      if( i > 0 ) os<<",";
      os<<s<<"="<<content[i].second;
    }
  os<<std::endl;
}

// Sets the value for the requested key
void CSparseSample::setValue( unsigned int key, double val )
{
  // Find the place for the value change
  std::vector< std::pair< unsigned int, double > >::iterator iter = content.begin();
  while( iter != content.end() && iter->first < key ) iter++;

  // Insertion
  if( iter == content.end() || iter->first > key )
    {
      if( val == 0.0 ) return;	// Maintain sparsity
      
      std::pair< unsigned int, double > entry( key, val );
      content.insert( iter, entry );
    }

  // Replacement
  else
    {
      if( val == 0.0 ) content.erase( iter );
      else iter->second = val;
    }
}

// Sets the value for the requested key
void CSparseSample::setValue( string key, double val )
{
  int j = pfmap->f2i( key );
  if( j < 0 ) throw std::logic_error( "No such feature" );
  setValue( j, val );
}

// Returns the value for the requested key
double CSparseSample::getValue( string key ) const
{
  int j = pfmap->f2i( key );
  if( j < 0 ) return 0.0;
  return getValue( j );
}

// Returns the value for the requested key
double CSparseSample::getValue( unsigned int key ) const
{
  for( unsigned int i = 0; i < content.size(); i++ )
    {
      if( content[i].first > key ) return 0.0;
      if( content[i].first == key ) return content[i].second;
    }
  return 0.0;
}

/// Maps the sample to a new feature map
void CSparseSample::remap( shared_ptr< CFeatMap const > newmap )
{
  // Generate the new content
  vector< std::pair< unsigned int, double > > new_content;
  for( unsigned int i = 0; i < newmap->nFeats(); i++ )
    {
      // Fetch the feature name
      string key = newmap->i2f( i );
      
      // Transfer the value
      double v = getValue( key );
      if( v == 0.0 ) continue;
      std::pair< unsigned int, double > entry( i, v );
      new_content.push_back( entry );
    }

  // Update the internal structures
  content.swap( new_content );
  pfmap = newmap;
}

// Computes the dot product with another sparse sample
double CSparseSample::dot( const CSparseSample& other ) const
{
  // Verify consistency of the feature maps
  if( pfmap != other.pfmap ) throw std::logic_error( "Incompatible computation in CSparseSample::dot()" );

  // Declarations
  double res = 0.0;
  unsigned int ii = 0;
  unsigned int ij = 0;
  unsigned int ni = this->content.size();
  unsigned int nj = other.content.size();
  
  // Traverse the samples
  while( ii < ni && ij < nj )
    {
      if( this->content[ii].first < other.content[ij].first ) ii++;
      else if( other.content[ij].first < this->content[ii].first ) ij++;
      else
	{
	  res += this->content[ii].second * other.content[ij].second;
	  ii++; ij++;
	}
    }

  return res;
}

// Multiplies the sparse vector by a scalar
CSparseSample& CSparseSample::operator*=( const double& coeff )
{
  // Degenerate case
  if( coeff == 0.0 ) {content.clear(); return *this;}
  
  // General case
  for( unsigned int i = 0; i < content.size(); ++i )
    { content[i].second *= coeff; }

  return *this;
}

// Adds another vector to *this
CSparseSample& CSparseSample::operator+=( const CSparseSample& rhs )
{
  std::vector< std::pair< unsigned int, double > > res;
  
  // Retrieve the size of each vector
  unsigned int ia = 0;
  unsigned int ib = 0;
  unsigned int na = this->content.size();
  unsigned int nb = rhs.content.size();

  // Traverse the vectors
  while( ia < na || ib < nb )
    {
      unsigned int key = 0;
      double value = 0.0;

      // Fetch the keys
      unsigned int keya = std::numeric_limits< unsigned int >::max();
      unsigned int keyb = std::numeric_limits< unsigned int >::max();

      if( ia < na ) keya = this->content[ia].first;
      if( ib < nb ) keyb = rhs.content[ib].first;

      // Fetch the value(s)
      if( keya < keyb )
	{
	  key = keya;
	  value = this->content[ia].second;
	  ++ia;
	}
      
      else if( keyb < keya )
	{
	  key = keyb;
	  value = rhs.content[ib].second;
	  ++ib;
	}

      else
	{
	  key = keya;	// == keyb
	  value = this->content[ia].second + 
	    rhs.content[ib].second;
	  ++ia;
	  ++ib;
	}

      // Store the value
      std::pair< unsigned int, double > p( key, value );
      res.push_back( p );
    }

  // By construction, res satisfies the class constraint
  // Store it
  this->content.swap( res );
  return *this;
}

// Multiples the sparse vector by a scalar
const CSparseSample CSparseSample::operator*( const double& coeff ) const
{
  CSparseSample res( *this );
  res *= coeff;
  return res;
}

// Key-based less-than
bool CSparseSample::operator<( const CSparseSample& other ) const
{
  // Retrieve the size of each vector
  unsigned int na = this->content.size();
  unsigned int nb = other.content.size();

  // Traverse the vectors
  for( unsigned int i = 0; i < na && i < nb; i++ )
    {
      std::pair< int, double > ea = this->content[i];
      std::pair< int, double > eb = other.content[i];

      if( ea.first < eb.first ) return true;
      if( eb.first < ea.first ) return false;
    }

  // Special case: sample a is shorter than pattern b
  if( na < nb ) return true;
  return false;
}

// Key-based equality
bool CSparseSample::operator==( const CSparseSample& other ) const
{
  // Retrieve the size of each vector
  unsigned int na = this->content.size();
  unsigned int nb = other.content.size();

  if( na != nb ) return false;

  // Traverse the vectors
  for( unsigned int i = 0; i < na && i < nb; i++ )
    {
      std::pair< int, double > ea = this->content[i];
      std::pair< int, double > eb = other.content[i];

      if( ea.first < eb.first ) return false;
      if( eb.first < ea.first ) return false;
    }
  
  return true;
}

// I/O
std::ostream& operator<<( std::ostream& os, const CSparseSample& smpl )
{
  smpl.display( os );
  return os;
}

// Returns the number of features in a dataset of sparse samples
unsigned int nFeats( const CDataSet<CSparseSample>& ds )
{
  if( ds.size() < 1 ) return 0;
  else return ds.getSample(0)->getFeatMap()->nFeats();
}

// Returns the number of features in a dataset of sparse samples
unsigned int nFeats( const shared_ptr< const CDataSet<CSparseSample> >& pds )
{
  if( !pds ) return 0;
  else return nFeats( *pds );
}

// Remaps the dataset to a new feature map
void remap( CDataSet<CSparseSample>& ds, shared_ptr< CFeatMap > fmap )
{
  // Project individual samples
  for( unsigned int i = 0; i < ds.size(); i++ )
    ds.getSampleMod(i)->remap(fmap);
}

// Removes all samples that have fewer non-zero features than thresh
void cropSamples( unsigned int thresh, CDataSet<CSparseSample>& ds )
{
  vector< unsigned int > toKeep;
  for( unsigned int i = 0; i < ds.size(); i++ )
    {
      unsigned int n = ds.getSample( i )->L0();
      if( n >= thresh ) toKeep.push_back( i );
    }
  ds.subsample( toKeep );
}

// Computes a set of common feature IDs between two datasets
vector< string > commonFeatIDs( const CDataSet<CSparseSample>& ds1,
				const CDataSet<CSparseSample>& ds2 )
{
  // Fetch the feature maps
  shared_ptr< CFeatMap const > pfm1 = ds1.getSample( 0 )->getFeatMap();
  shared_ptr< CFeatMap const > pfm2 = ds2.getSample( 0 )->getFeatMap();

  // Fetch the ids
  vector< string > v1 = pfm1->getFeatureIDs();
  vector< string > v2 = pfm2->getFeatureIDs();

  // Compute the intersection
  vector< string > res;
  std::sort( v1.begin(), v1.end() );
  std::sort( v2.begin(), v2.end() );
  std::set_intersection( v1.begin(), v1.end(), v2.begin(), v2.end(),
			 std::inserter( res, res.begin() ) );
  return res;
}

// Returns the range of values in the dataset
std::pair< double, double > getRange( const CDataSet<CSparseSample>& ds )
{
  double min = std::numeric_limits<double>::infinity();
  double max = -1.0 * min;
  for( unsigned int i = 0; i < ds.size(); i++ )
    {
      shared_ptr< CSparseSample const > s = ds.getSample( i );
      for( unsigned int j = 0; j < nFeats(ds); j++ )
	{
	  double v = s->getValue( j );
	  if( v > max ) max = v;
	  else if( v < min ) min = v;
	}
    }
  return std::pair<double, double>( min, max );
}

// Computes the average precision and recall per feature, given the truth
std::pair< double, double > computePnR( const CDataSet<CSparseSample>& dsPred, const CDataSet<CSparseSample>& dsTruth )
{
  // Sums of precisions and recalls
  double psum = 0.0;
  double rsum = 0.0;

  // Retrieve the feature maps
  shared_ptr< CFeatMap const > pfmTruth = dsTruth.getSample(0)->getFeatMap();
  shared_ptr< CFeatMap const > pfmPred = dsPred.getSample(0)->getFeatMap();

  // Traverse the features
  for( unsigned int iFeat = 0; iFeat < nFeats(dsPred); iFeat++ )
    {
      unsigned int nP = 0;
      unsigned int nTP = 0;
      unsigned int nFP = 0;

      // Find the feature in the "truth"
      int jFeat = pfmTruth->f2i( pfmPred->i2f(iFeat) );
      if( jFeat < 0.0 ) throw std::logic_error( "computePnR(): No such feature" );

      for( unsigned int i = 0; i < dsPred.size(); i++ )
	{
	  // Find the sample in the "truth"
	  int j = dsTruth.s2i( dsPred.i2s(i) );
	  if( j < 0.0 ) throw std::logic_error( "computePnR(): No such sample" );
	  
	  // Retrieve the prediction and the truth
	  bool bTruth = dsTruth.getSample(j)->getValue( jFeat ) > 0.0;
	  bool bPred = dsPred.getSample(i)->getValue( iFeat ) > 0.0;

	  if( bTruth ) nP++;
	  if( bPred )
	    {
	      if( bTruth ) nTP++;
	      else nFP++;
	    }
	}

      // Compute precision and recall
      double p = 0.0;
      if( nTP + nFP > 0 )
	p = static_cast<double>( nTP ) / static_cast<double>( nTP + nFP );
      double r = static_cast<double>( nTP ) / static_cast<double>( nP );

      psum += p;
      rsum += r;
    }

  // Compute and return the averages
  std::pair< double, double > res;
  res.first = psum / static_cast<double>( nFeats(dsPred) );
  res.second = rsum / static_cast<double>( nFeats(dsPred) );
  return res;
}

// Thresholds the dataset, setting the top-scoreing k samples for each feature to 1 and rest to 0
void threshold( CDataSet<CSparseSample>& ds, std::map< string, unsigned int > prof )
{
  // Fetch the feature map
  shared_ptr< CFeatMap const > pfm = ds.getSample(0)->getFeatMap();

  // Traverse the features
  for( unsigned int iFeat = 0; iFeat < nFeats(ds); iFeat++ )
    {
      // Find the feature in the provided profile
      string s = pfm->i2f( iFeat );
      std::map< string, unsigned int >::iterator iter = prof.find( s );
      if( iter == prof.end() ) throw std::logic_error( "threshDataSet: Invalid profile" );
      unsigned int k = iter->second;

      // Set-up alphabetical access to the samples
      // This ensures invariance to subsample(random permutation)
      vector< string > sids;
      for( unsigned int i = 0; i < ds.size(); i++ )
	sids.push_back( ds.i2s(i) );
      std::sort( sids.begin(), sids.end() );

      // Rank the values using alphabetical access
      std::multimap< double, unsigned int > m;
      for( unsigned int i = 0; i < sids.size(); i++ )
	{
	  std::pair<double, unsigned int> entry;
	  int j = ds.s2i(sids[i]);
	  entry.first = ds.getSample(j)->getValue( iFeat );
	  entry.second = j;
	  m.insert( entry );
	}

      // Traverse the values from highest to lowest
      unsigned int nPos = 0;
      for( std::multimap< double, unsigned int >::reverse_iterator iter = m.rbegin();
	   iter != m.rend(); iter++ )
	{
	  shared_ptr< CSparseSample > s = ds.getSampleMod( iter->second );
	  if( nPos < k )
	    {
	      s->setValue( iFeat, 1.0 );
	      nPos++;
	    }
	  else
	    s->setValue( iFeat, 0.0 );
	}
    }
}

// Thresholds the dataset, setting all entries for feature k that exceed prof[k] to 1 and rest to 0
void threshold( CDataSet<CSparseSample>& ds, std::map< string, double > prof )
{
  // Fetch the feature map
  shared_ptr< CFeatMap const > pfm = ds.getSample( 0 )->getFeatMap();

  // Traverse the features
  for( unsigned int iFeat = 0; iFeat < nFeats(ds); iFeat++ )
    {
      // Find the feature in the provided profile
      string s = pfm->i2f( iFeat );
      std::map< string, double >::iterator iter = prof.find( s );
      if( iter == prof.end() ) throw std::logic_error( "threshDataSet: Invalid profile" );
      double thr = iter->second;

      // Traverse the values
      for( unsigned int i = 0; i < ds.size(); i++ )
	{
	  shared_ptr< CSparseSample > s = ds.getSampleMod( i );
	  if( s->getValue( iFeat ) > thr )
	    s->setValue( iFeat, 1.0 );
	  else
	    s->setValue( iFeat, 0.0 );
	}
    }
}

// Returns the number of samples that have feature i
unsigned int nSamplesWFeat( const CDataSet<CSparseSample>& ds, string f )
{
  unsigned int n = 0;
  for( unsigned int j = 0; j < ds.size(); j++ )
    {
      if( ds.getSample(j)->getValue( f ) != 0.0 ) n++;
    }
  return n;
}

/// Returns a mapping of features to the number of samples that have it
void nSamplesWFeat( const CDataSet<CSparseSample>& ds, shared_ptr< CFeatMap const > pfmap, simap_t& counts )
{
  for( unsigned int j = 0; j < pfmap->nFeats(); j++ )
    {
      // Retrieve the feature name
      string f = pfmap->i2f( j );
      
      // Count the number of samples with this feature in the dataset
      unsigned int n = nSamplesWFeat( ds, f );
	    
      // Append the count to counts
      simap_t::iterator iter = counts.find( f );
      if( iter == counts.end() ) counts[f] = n;
      else counts[f] += n;
    }
}

// Adds new features to the sameples of ds, overwriting any previous values
void expand( CDataSet<CSparseSample>& ds, shared_ptr< CFeatMap > pfm,
	     const CDataSet<CSparseSample>& other )
{
  // Retrieve the new feature map
  shared_ptr< const CFeatMap > pfmOther = other.getSample(0)->getFeatMap();

  // Traverse the new features and add them to pfm
  for( unsigned int i = 0; i < pfmOther->nFeats(); i++ )
    {
      string str = pfmOther->i2f( i );
      pfm->addFeat( str );
    }

  // Add individual feature values to the dataset
  for( unsigned int i = 0; i < other.size(); i++ )
    {
      // Locate the sample in ds
      string name = other.i2s( i );
      int iSmpl = ds.s2i( name );

      // If it doesn't exist, create it
      if( iSmpl < 0 )
	{
	  auto_ptr< CSparseSample > ps( new CSparseSample( pfm ) );
	  ds.addSample( name, ps );
	  iSmpl = ds.s2i( name );
	}
      if( iSmpl < 0 ) throw std::logic_error( "Failed to locate a sample after adding it" );
      
      // Traverse the features
      for( unsigned int j = 0; j < pfmOther->nFeats(); j++ )
	{
	  // Retrieve the feature name and its value
	  string fName = pfmOther->i2f( j );
	  double fVal = other.getSample( i )->getValue( j );

	  // Update ds
	  int jFeat = pfm->f2i( fName );
	  ds.getSampleMod( iSmpl )->setValue( jFeat, fVal );
	}
    }
}

// Adds new features to the samples of ds, appending them as another feature space
void expand( CDataSet<vSparseSample>& ds, const CDataSet<CSparseSample>& other, bool bRemoveMissing )
{
  // Handle the degenerate case
  if( other.size() == 0 ) return;

  // Retrieve the current number of kernel spaces in ds
  unsigned int nSpaces = 0;
  if( ds.size() > 0 ) nSpaces = ds.getSample( 0 )->size();

  // Traverse the new samples
  vector< string > toKeep;
  for( unsigned int i = 0; i < other.size(); ++i )
    {
      // Locate the sample in ds
      string name = other.i2s( i );
      int iSmpl = ds.s2i( name );

      // Handle the missing data
      if( iSmpl < 0 )
	{
	  if( bRemoveMissing ) continue;

	  // Add the empty samples for the current kernel spaces
	  auto_ptr< vSparseSample > pv( new vSparseSample );
	  for( unsigned int j = 0; j < nSpaces; ++j )
	    {
	      // nSpaces > 0 guarantees that there's at least one example in ds; retrieve it
	      shared_ptr< const vSparseSample > pv0 = ds.getSample( 0 );

	      // Retrieve the feature map associated with the kernel space j
	      shared_ptr< const CFeatMap > pfm = (*pv0)[j].getFeatMap();

	      // Add an empty sample for the kernel space j
	      CSparseSample smpl( pfm );
	      pv->push_back( smpl );
	    }
	  ds.addSample( name, pv );
	  iSmpl = ds.s2i( name );
	}
      if( iSmpl < 0 ) throw std::logic_error( "Failed to locate a sample after adding it" );

      // Append the new feature space
      CSparseSample smpl( *other.getSample(i) );
      ds.getSampleMod( iSmpl )->push_back( smpl );
      toKeep.push_back( name );
    }

  // Handle the intersection
  if( bRemoveMissing ) ds.subsample( toKeep );

  // Handle the union
  else
    {
      // Pad the missing data with empty samples
      shared_ptr< const CFeatMap > pfmOther = other.getSample( 0 )->getFeatMap();
      for( unsigned int i = 0; i < ds.size(); ++i )
	{
	  shared_ptr< vSparseSample > pv = ds.getSampleMod( i );
	  if( pv->size() < nSpaces+1 )
	    {
	      CSparseSample smpl( pfmOther );
	      pv->push_back( smpl );
	    }
	}
    }
}

// Returns the number of features in an MKL dataset
unsigned int nFeats( const CDataSet<vSparseSample>& ds )
{
  if( ds.size() < 1 ) return 0;
  unsigned int res = 0;
  shared_ptr< const vSparseSample > ps = ds.getSample( 0 );
  for( unsigned int i = 0; i < ps->size(); ++i )
    res += (*ps)[i].getFeatMap()->nFeats();
  return res;
}

/// Returns the number of kernels used by an MKL dataset
unsigned int nKernels( const CDataSet<vSparseSample>& ds )
{
  if( ds.size() < 1 ) return 0;
  else return ds.getSample(0)->size();
}

// Takes a set of protein names and creates a sparse dataset over the annotations for those proteins
void makeSparseDataset( const GO::CGOACollection& goaSource,
			const vector<string>& protnames,
			CDataSet<CSparseSample>& ds,
			const GO::CGOContainer& goGraph,
			GO::ONTOLOGY_INDEX filter,
			shared_ptr< CFeatMap > pfmap )
{
  // Create a feature map if none was given
  if( !pfmap )
    pfmap.reset( new CFeatMap );

  // Traverse the requested proteins
  for( vector< string >::const_iterator iter = protnames.begin();
       iter != protnames.end(); iter++ )
    {
      // Retrieve the set of GO IDs for the protein and project them onto the slim ontology
      set< string > annot_full;
      vector< string > annots = goaSource.getGOIDs( *iter, filter );

      // Don't add proteins with no annotations
      if( annots.size() < 1 ) continue;

      goGraph.getFullPaths( annots, annot_full );

      // Traverse the annotations and compose the sample
      std::auto_ptr< CSparseSample > ps( new CSparseSample( pfmap ) );
      for( std::set< std::string >::iterator annot_iter = annot_full.begin();
	   annot_iter != annot_full.end(); annot_iter++ )
	{
	  // Retrieve the feature index
	  unsigned int fi = pfmap->addFeat( *annot_iter );

	  // Add the key/value pair to the sample
	  ps->setValue( fi, 1.0 );
	}

      // Add the sample
      ds.addSample( *iter, ps );
    }
}

void makeSparseDataset( const GO::CBLASTOutput& source,
			CDataSet<CSparseSample>& ds,
			double lower_thresh, double upper_thresh,
			shared_ptr< CFeatMap > pfmap )
{
  typedef GO::CBLASTOutput::SBOEntry SBOEntry;
  typedef GO::CBLASTOutput::emap_const_iter_t iter_t;

  // Create a feature map if none was given
  if( !pfmap )
    pfmap.reset( new CFeatMap );

  // Traverse the requested proteins
  for( iter_t iter = source.begin(); iter != source.end(); ++iter )
    {
      // Retrieve the hits
      const vector< SBOEntry >& vHits = iter->second;
      map< string, double > mHits;
      for( unsigned int k = 0; k < vHits.size(); k++ )
	{
	  // Retrieve the key/value pair
	  string key = vHits[k].subject_id;
	  double val = vHits[k].e_value;
	  //	  double pid = vHits[k].percent_identity;

	  // Skip self-hit
	  if( key == iter->first ) continue;

	  // Skip 100% identity hits
	  //	  if( pid == 100.0 ) continue;
	  
	  // Skip large e-values
	  if( val > upper_thresh ) continue;

	  // Determine if the key already exists in the map
	  if( mHits.find( key ) != mHits.end() )
	    {
	      // If the new e-value is greater than the previous one, do nothing
	      if( mHits[key] < val ) continue;
	    }
	      
	  // Otherwise, insert/overwrite it
	  mHits[key] = val;
	}
      
      // Do nothing if there are no hits that satisfy the filter
      if( mHits.empty() ) continue;

      // The map ensures that there's only one entry for each key
      // Compose the sample
      auto_ptr< CSparseSample > ps( new CSparseSample( pfmap ) );

      // Traverse the map and insert the features
      for( map< string, double >::iterator hiter = mHits.begin();
	   hiter != mHits.end(); ++hiter )
	{
	  // Find the feature in the map
	  unsigned int fi = pfmap->addFeat( hiter->first );
	
	  // Compute the log-score
	  double d = hiter->second / upper_thresh;
	  if( d < lower_thresh ) d = lower_thresh;
	  d = -1.0 * log( d );

	  // Add the key/value pair to the sample
	  ps->setValue( fi, d );
	}

      // Add the entry to the dataset
      ds.addSample( iter->first, ps );
    }
}
