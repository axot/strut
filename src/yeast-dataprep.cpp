// yeast-dataprep.cpp - Prepares data for yeast experiments
//
// by Artem Sokolov

#include "types.h"
#include "parsers.h"
#include "blastout.h"
#include "go-annotation.h"

#include <string>
#include <iostream>
#include <fstream>

#include<boost/tokenizer.hpp>

using std::cout;
using std::endl;
using std::ifstream;

const GO::ONTOLOGY_INDEX g_myfilter = GO::GO_MF;
const string fnPrefix = "yeast/mf_";

const string blast_filename( "/s/chopin/c/proj/protfun/data/BLAST/s_cerevisiae/s_cerevisiae-2009-06-04.blast" );
const string ppi_filename( "/s/chopin/c/proj/protfun/data/organisms/s_cerevisiae/ppi/biogrid.data" );
const string tm_filename( "/s/chopin/c/proj/protfun/data/organisms/s_cerevisiae/yeast_tm.data" );
const string loc_filename( "/s/chopin/c/proj/protfun/data/organisms/s_cerevisiae/yeast_localization.data" );
const string annot_filename( "/s/chopin/c/proj/protfun/data/organisms/s_cerevisiae/s_cerevisiae-2009-06-04.annot.gz" );

// Set intersection between two sorted vectors with an intuitive interface
template< typename _T >
vector<_T> set_intersection( const vector<_T>& u, const vector<_T>& v )
{
  vector<_T> res;
  std::set_intersection( u.begin(), u.end(), v.begin(), v.end(),
			 std::inserter( res, res.begin() ) );
  return res;
}

// Renames the sample ids of a dataset based on the provided annotation collection
// If the id is not in the collection, it is not modified
// Returns the set of modified ids
const vector< string > renameIDs( shared_ptr< CDataSet<CSparseSample> > pds, const GO::CGOACollection& goa )
{
  vector< string > res;

  // Traverse the samples
  for( unsigned int i = 0; i < pds->size(); i++ )
    {
      string s = pds->i2s( i );
      string id = goa.getObjectID( s );
      if( id.size() > 0 )
	{
	  try { pds->rename( i, id ); res.push_back( id ); }
	  catch( std::logic_error e )
	    {
	      cout<<"Dealing with sample "<<i<<endl;
	      cout<<"Trying to rename "<<s<<" to "<<id<<endl;
	      cout<<"Name "<<id<<" is already in the dataset as sample "<<pds->s2i( id )<<endl;
	    }
	}
    }
  
  cout<<"Renamed "<<res.size()<<" samples"<<endl;
  return res;
}

shared_ptr< CDataSet<CSparseSample> > parseLocalization( string filename )
{
  typedef boost::tokenizer< boost::char_separator<char> > tokenizer;
  boost::char_separator<char> spr( "," );

  // Open the file
  ifstream ifs( filename.c_str() );
  if( ifs.fail() )
    throw std::runtime_error( "Failed to open " + filename );

  shared_ptr< CDataSet<CSparseSample> > pRes( new CDataSet<CSparseSample> );
  string s;

  // Parse the header
  vector< string > featIDs;
  std::getline( ifs, s );
  tokenizer tk( s, spr );
  tokenizer::iterator iter = tk.begin();
  ++iter; ++iter;	// Skip "patternID" and "labels"
  for( ; iter != tk.end(); ++iter )
    featIDs.push_back( *iter );

  cout<<featIDs.size()<<" features... "; cout.flush();

  // Create the feature map
  shared_ptr< CFeatMap > pfm( new CFeatMap( featIDs ) );

  // Parse the samples
  std::getline( ifs, s );
  while( ifs.fail() == false )
    {
      if( s.size() == 0 ) break;

      string id;
      vector< double > vals;

      // Parse the sampleID and the values
      tokenizer tk( s, spr );
      tokenizer::iterator iter = tk.begin();
      id = *iter; ++iter;	// sample ID
      ++iter;			// skip 'unknown'
      for( ; iter != tk.end(); ++iter )
	vals.push_back( boost::lexical_cast<double>( *iter ) );

      // Verify the dimensionality
      if( vals.size() != featIDs.size() )
	throw std::logic_error( "Number of values doesn't match the number of features" );

      // Create the sample
      auto_ptr< CSparseSample > ps( new CSparseSample(pfm) );
      for( unsigned int i = 0; i < vals.size(); i++ )
	ps->setValue( featIDs[i], vals[i] );

      // Add the sample
      pRes->addSample( id, ps );

      std::getline( ifs, s );
    }

  return pRes;
}

int main( int argc, char* argv[] )
{
  if( argc < 2 )
    {
      cout<<"Usage: "<<argv[0]<<" <location of the .obo file> "<<endl;
      return -1;
    }

  // Preload the ontology
  string obo_loc( argv[1] );
  GO::CGOContainer goGraph( obo_loc );

  cout<<"Loading the BLAST hits data... "; cout.flush();
  GO::CBLASTOutput blastout( blast_filename.c_str() );
  cout<<blastout.size()<<" entires loaded"<<endl;

  cout<<"Loading annotations... "; cout.flush();
  GO::CGOACollection goa( annot_filename.c_str() );
  cout<<goa.size()<<" annotations loaded"<<endl;

  cout<<"Generating a BLAST dataset... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsBlast( new CDataSet<CSparseSample> );
  makeSparseDataset( blastout, *pdsBlast, 1e-10, 50.0 );
  cout<<pdsBlast->size()<<" samples generated"<<endl;

  cout<<"Loading transmembrane data... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsTM( new CDataSet<CSparseSample> );
  parseSparseFile( tm_filename, *pdsTM, ',', ' ', ':' );
  cout<<pdsTM->size()<<" samples loaded"<<endl;

  cout<<"Loading localization data... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsLoc = parseLocalization(loc_filename);
  cout<<pdsLoc->size()<<" samples loaded"<<endl;

  cout<<"Loading PPI data... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsPPI( new CDataSet<CSparseSample> );
  parseSparseFile( ppi_filename, *pdsPPI, ',', ' ', ':' );
  cout<<pdsPPI->size()<<" samples loaded"<<endl;

  // Determine the intersection of sample IDs
  vector< string > v1 = pdsBlast->getSampleIDs();
  vector< string > v2 = pdsTM->getSampleIDs();
  vector< string > v3 = pdsLoc->getSampleIDs();
  vector< string > v4 = renameIDs( pdsPPI, goa );

  std::sort( v1.begin(), v1.end() );
  std::sort( v2.begin(), v2.end() );
  std::sort( v3.begin(), v3.end() );
  std::sort( v4.begin(), v4.end() );

  // Determine the common set of ids
  vector< string > v12 = set_intersection( v1, v2 );
  cout<<v12.size()<<" sample IDs are common to 1 and 2"<<endl;
  vector< string > v123 = set_intersection( v12, v3 );
  cout<<v123.size()<<" sample IDs are common to 1, 2 and 3"<<endl;
  vector< string > v1234 = set_intersection( v123, v4 );
  cout<<v1234.size()<<" sample IDs are common to all datasets"<<endl;

  cout<<"Generating the annotation dataset... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsAnnot( new CDataSet<CSparseSample> );
  makeSparseDataset( goa, v1234, *(pdsAnnot), goGraph, g_myfilter );
  cout<<pdsAnnot->size()<<" samples; "<<nFeats( pdsAnnot )<<" features generated"<<endl;
  shared_ptr< CFeatMap const > pfm = pdsAnnot->getSample(0)->getFeatMap();

  // Isolate well-represented features
  const int reprThreshold = 10;
  simap_t featCount;
  nSamplesWFeat( *(pdsAnnot), pfm, featCount );
  vector< string > reprFeat;
  for( simap_t::const_iterator iter = featCount.begin(); iter != featCount.end(); iter++ )
    if( iter->second >= reprThreshold ) reprFeat.push_back( iter->first );
  shared_ptr< CFeatMap > pfmRepr( new CFeatMap( reprFeat ) );
  remap( *pdsAnnot, pfmRepr );

  // Isolate well-represented samples
  cropSamples( 2, *pdsAnnot );
  cout<<pdsAnnot->size()<<" samples; "<<nFeats( pdsAnnot )<<" features are well-represented"<<endl;
 
  // Reduce each dataset to the common set of sample IDs
  vector< string > vids = pdsAnnot->getSampleIDs();
  pdsBlast->subsample( vids );
  pdsTM->subsample( vids );
  pdsLoc->subsample( vids );
  pdsPPI->subsample( vids );

  // Save the datasets
  string fnBlast = fnPrefix + "blast.sdat"; pdsBlast->save( fnBlast );
  string fnTM = fnPrefix + "transmem.sdat"; pdsTM->save( fnTM );
  string fnLoc = fnPrefix + "local.sdat"; pdsLoc->save( fnLoc );
  string fnPPI = fnPrefix + "ppi.sdat"; pdsPPI->save( fnPPI );
  string fnAnnot = fnPrefix + "annot.sdat"; pdsAnnot->save( fnAnnot );
  
  return 0;
}
