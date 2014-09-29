// cafa-dataprep.cpp - preparation of datasets for the CAFA challenge
//
// by Artem Sokolov

#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "misc.h"
#include "types.h"
#include "go-annotation.h"
#include "blastout.h"
#include "sample.h"
#include "parsers.h"

namespace po = boost::program_options;

typedef struct
{
  string name;
  double pIdent;
  unsigned int mLength;
} id_entry;

typedef std::unordered_map< string, id_entry > idmap_t;
typedef std::unordered_map< string, string > ssmap_t;

// Set intersection intuitive interface
template< typename _T >
set<_T> set_intersection( const set<_T>& u, const set<_T>& v )
{
  set<_T> res;
  std::set_intersection( u.begin(), u.end(), v.begin(), v.end(),
			 std::inserter( res, res.begin() ) );
  return res;
}

// Set difference intuitive interface
template< typename _T >
set<_T> set_difference( const set<_T>& u, const set<_T>& v )
{
  set<_T> res;
  std::set_difference( u.begin(), u.end(), v.begin(), v.end(),
		       std::inserter( res, res.begin() ) );
  return res;
}

// Appends dataset B to dataset A
void append( CDataSet<CSparseSample>& A, const CDataSet<CSparseSample>& B, shared_ptr< CFeatMap > pfmA )
{
  // Degenerate case
  if( B.size() == 0 ) return;

  // Fetch the feature map of B
  shared_ptr< CFeatMap const > pfmB = B.getSample( 0 )->getFeatMap();

  // Traverse the samples
  for( unsigned int i = 0; i < B.size(); ++i )
    {
      string s = B.i2s( i );

      // Source
      shared_ptr< CSparseSample const > ps = B.getSample( i );

      // Destination
      int j = A.s2i( s );
      if( j < 0 )
	{
	  auto_ptr< CSparseSample > p( new CSparseSample( pfmA ) );
	  A.addSample( s, p );
	  j = A.s2i( s );
	  if( j < 0 ) throw std::logic_error( "Sample does not exist after its addition to the dataset" );
	}
      shared_ptr< CSparseSample > pd = A.getSampleMod( static_cast<unsigned int>(j) );
      
      // Traverse the features
      for( unsigned int fi = 0; fi < pfmB->nFeats(); ++fi )
	{
	  string f = pfmB->i2f( fi );
	  double val = ps->getValue( f );
	  if( val != 0.0 )
	    {
	      pfmA->addFeat( f );
	      pd->setValue( f, val );
	    }
	}
    }
}

// Loads a dataset, merges it with the joint one, and displays some statistics
void load_and_append( const string& name, const string& filename, CDataSet<CSparseSample>& dsJoint, shared_ptr< CFeatMap > pfm, char sep1 = ' ', char sep2 = ':' )
{
  CDataSet<CSparseSample> ds;

  // Load the data
  parseSparseFile( filename, ds, ',', sep1, sep2 );
  if( ds.size() == 0 || nFeats( ds ) == 0 ) throw std::logic_error( "Failed to load " + name + " data" );
  cout<<"  "<<name<<" data has "<<ds.size()<<" samples; ";
  cout<<nFeats( ds )<<" features"<<endl;
  cout<<"    Example id: "<<ds.i2s(0)<<endl;

  // Append the data to the joint dataset
  append( dsJoint, ds, pfm );
  cout<<"  Joint data has "<<dsJoint.size()<<" samples; ";
  cout<<nFeats( dsJoint )<<" features"<<endl;
}

// Loads an input space representation
CDataSet<CSparseSample> load_ispace( const string& fnPrefix, shared_ptr< CFeatMap > pfm )
{
  CDataSet< CSparseSample > dsRes;
  CDataSet< CSparseSample > dsB;
  CDataSet< CSparseSample > dsLC;
  CDataSet< CSparseSample > dsLZ;
  CDataSet< CSparseSample > dsTN;
  CDataSet< CSparseSample > dsTM;

  // Compose filenames using the prefix
  string fnB = fnPrefix + ".blast.gz";
  string fnLC = fnPrefix + ".lcomplx.dat";
  string fnLZ = fnPrefix + ".localiz.dat";
  string fnTN = fnPrefix + ".termini.dat";
  string fnTM = fnPrefix + ".transmem.dat";

  // Load the BLAST data
  cout<<"  Loading the BLAST training data... "; cout.flush();
  GO::CBLASTOutput bo( fnB.c_str() );
  cout<<"loaded "<<bo.size()<<" entries"<<endl;

  // Generate the BLAST dataset
  makeSparseDataset( bo, dsRes, 1e-10, 50.0, pfm );
  if( dsRes.size() == 0 || nFeats( dsRes ) == 0 ) throw std::logic_error( "Failed to load the BLAST data" );
  cout<<"  BLAST data has "<<dsRes.size()<<" samples; ";
  cout<<nFeats( dsRes )<<" features"<<endl;
  cout<<"    Example id: "<<dsRes.i2s(0)<<endl;

  // Load all other dataset
  load_and_append( "low-complexity", fnLC, dsRes, pfm );
  load_and_append( "localization", fnLZ, dsRes, pfm );
  load_and_append( "termini", fnTN, dsRes, pfm );
  load_and_append( "transmem", fnTM, dsRes, pfm );

  return dsRes;
}

// Loads the ID map
idmap_t loadIDMap( const string& filename )
{
  idmap_t res;

  // Open the file and traverse it line-by-line
  std::ifstream ifs( filename.c_str() );
  string s;
  for( std::getline( ifs, s ); ifs.fail() == false; std::getline( ifs, s ) )
    {
      // Tokenize the string
      vector< string > toks;
      boost::split( toks, s, boost::is_space() );
      if( toks.size() < 4 ) throw std::logic_error( "Invalid format of the id match file" );

      // Parse percent-identy and match-length
      double pi = boost::lexical_cast<double>( toks[2] );
      double ml = boost::lexical_cast<unsigned int>( toks[3] );

      // Handle existing entry
      if( res.find( toks[0] ) != res.end() )
	{
	  // Percent identity needs to be higher for replacement
	  if( pi < res[toks[0]].pIdent ) continue;

	  // Match length needs to be higher for replacement
	  if( ml < res[toks[0]].mLength ) continue;
	}

      // Add / replace the entry
      id_entry e = { toks[1], pi, ml };
      res[toks[0]] = e;
    }

  return res;
}

// Inverts the id map
ssmap_t invertIDMap( const idmap_t& m )
{
  ssmap_t res;
  for( idmap_t::const_iterator iter = m.begin(); iter != m.end(); ++iter )
    {
      string key = iter->second.name;
      string val = iter->first;
      res[key] = val;
    }
  return res;
}

// Maps a set of IDs in a dataset using the provided ID map
set< string > mapIDs( const CDataSet<CSparseSample>& ds, const string& fn )
{
  // Load the id map
  cout<<"Loading the id map... "; cout.flush();
  idmap_t m = loadIDMap( fn );
  cout<<"loaded "<<m.size()<<" entries"<<endl;

  // Traverse the dataset IDs and use the mapping
  set< string > res;
  for( unsigned int i = 0; i < ds.size(); ++i )
    {
      string s = ds.i2s( i );
      idmap_t::const_iterator iter = m.find( s );
      if( iter == m.end() ) continue;
      res.insert( iter->second.name );
    }
  return res;
}

// Removes any examples with names in rm and returns a set of common ids to an input space sampling
vector< string > filter_ispace( CDataSet<CSparseSample>& is, const set<string>& rm )
{
  // Find the set of common ids
  vector< string > vInit = is.getSampleIDs();
  set< string > sCommon( vInit.begin(), vInit.end() );
  cout<<sCommon.size()<<" ids are in the training data"<<endl;
  if( sCommon.empty() == false )
    cout<<"Example of a common id: "<<*sCommon.begin()<<endl;
  
  // Remove the prohibited ids
  set< string > sTemp = set_difference( sCommon, rm );
  vector< string > v( sTemp.begin(), sTemp.end() );

  // Subsample the datasets according to the resulting set of ids
  is.subsample( v );

  // Sanity checks
  for( vector< string >::iterator iter = v.begin(); iter != v.end(); ++iter )
    {
      if( is.s2i( *iter ) < 0 )
	throw std::logic_error( "Invalid filtering" );
    }
  for( unsigned int i = 0; i < is.size(); ++i )
    {
      string s = is.i2s( i );
      if( rm.find( s ) != rm.end() )
	throw std::logic_error( "Prohibited ids still exist" );
    }

  return v;
}

// Handles PPI data
void process_ppi( const string& fnPPI, const set<string>& toRemove, const string& sOPrefix, const string& fnID )
{
  cout<<"Processing PPI data"<<endl;

  cout<<"  Loading the id map for PPI data processing... "; cout.flush();
  idmap_t mOrig = loadIDMap( fnID );
  cout<<"loaded "<<mOrig.size()<<" entries"<<endl;

  // Invert the ID map
  ssmap_t m = invertIDMap( mOrig );
  
  // Sanity check
  for( ssmap_t::iterator iter = m.begin(); iter != m.end(); ++iter )
    {
      string key = iter->second;
      string val = iter->first;
      if( mOrig[key].name != val )
	throw std::logic_error( "Invalid ID map inversion" );
    }

  cout<<"  Loading PPI data... "; cout.flush();
  shared_ptr< CFeatMap > pfm( new CFeatMap );
  CDataSet<CSparseSample> dsPPI;
  parseSparseFile( fnPPI, dsPPI, pfm, ',', ',', '=' );
  cout<<"parsed "<<dsPPI.size()<<" samples, "<<nFeats( dsPPI )<<" features"<<endl;

  // Special case: Mouse data
  if( fnPPI.find( "Euk4" ) != fnPPI.npos )
    {
      string fnCOM = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/euk/Euk4/co-mention/co-mention.sdat.gz";
      load_and_append( "co-mention", fnCOM, dsPPI, pfm, ',', '=' );
    }

  // Retrieve the sample IDs and separate them into training and test sets
  vector< string > vAll = dsPPI.getSampleIDs();
  set< string > sAll( vAll.begin(), vAll.end() );
  set< string > sTrain = set_difference( sAll, toRemove );
  vector< string > vTrain( sTrain.begin(), sTrain.end() );
  vector< string > vTest( toRemove.begin(), toRemove.end() );  
  
  // Open the files
  boost::iostreams::filtering_ostream ofsPPITr;
  boost::iostreams::filtering_ostream ofsPPITs;
  openWriteFile( sOPrefix + "ppi_train.sdat.gz", ofsPPITr );
  openWriteFile( sOPrefix + "ppi_test.sdat.gz", ofsPPITs );

  // Save the training part
  dsPPI.displaySamples( vTrain, ofsPPITr );

  // Save the test samples using the original test IDs
  //  dsPPI.displaySamples( vTest, ofsPPITs );
  for( unsigned int i = 0; i < vTest.size(); ++i )
    {
      // Fetch the name
      string s = m[ vTest[i] ];

      // Fetch the sample
      int j = dsPPI.s2i( vTest[i] );
      if( j < 0 ) continue;
      shared_ptr< const CSparseSample > ps = dsPPI.getSample( j );

      // Write the sample
      ofsPPITs<<s<<","<<*ps;
    }
}

// Entry point
int main( int argc, char* argv[] )
{
  string fnIDs( "" );
  string fnOBO( "" );
  string fnAnnots( "" );
  string fnPPI( "" );
  string sOPrefix( "" );
  string sTrPrefix( "" );
  string sTsPrefix( "" );

  // Compose the allowed program arguments
  po::options_description opts( "Supported options" );
  opts.add_options()
    ( "obo-location,l", po::value< string >( &fnOBO ), "Ontology file to be used with goa files" )
    ( "id-match,m", po::value< string >( &fnIDs ), "Location of the file matching test IDs" )
    ( "annots,n", po::value< string >( &fnAnnots ), "Location of the annotations file" )
    ( "output-prefix,o", po::value< string >( &sOPrefix ), "Prefix for output files" ) 
    ( "ppi-data,p", po::value< string >( &fnPPI ), "File containing the PPI data" )
    ( "test-prefix,s", po::value< string >( &sTsPrefix ), "Prefix for test input files" )
    ( "train-prefix,t", po::value< string >( &sTrPrefix ), "Prefix for training input files" )
    ;

  // Display usage as needed
  if( argc < 2 )
    {
      cout<<"Usage: "<<argv[0]<<" options"<<endl;
      cout<<opts;
      return -1;
    }

  // Parse the arguments
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, opts), vm);
  po::notify(vm);

  // Display the parsed arguments
  cout<<".obo file: "<<fnOBO<<endl;
  cout<<"ID match file: "<<fnIDs<<endl;
  cout<<"Annotations: "<<fnAnnots<<endl;
  cout<<"PPI data: "<<fnPPI<<endl;
  cout<<"Output prefix: "<<sOPrefix<<endl;
  cout<<"Training data prefix: "<<sTrPrefix<<endl;
  cout<<"Test data prefix: "<<sTsPrefix<<endl;

  // Common feature maps
  shared_ptr< CFeatMap > pfm( new CFeatMap );

  // Load the training data
  cout<<"Loading training data"<<endl;
  CDataSet<CSparseSample> isTrain;
  if( sTrPrefix.empty() == false ) isTrain = load_ispace( sTrPrefix, pfm );
  vector< string > vGood = isTrain.getSampleIDs();
  set< string > toRemove;

  // Handle the test data
  if( sTsPrefix.empty() == false )
    {
      cout<<"Loading test data"<<endl;
      CDataSet<CSparseSample> isTest = load_ispace( sTsPrefix, pfm );

      // Save the data
      boost::iostreams::filtering_ostream ofsTs;
      openWriteFile( sOPrefix + "test.sdat.gz", ofsTs );
      isTest.display( ofsTs );

      // Find a set of ids matched to the test data
      if( fnIDs.empty() == false )
	{
	  toRemove = mapIDs( isTest, fnIDs );
	  cout<<toRemove.size()<<" test proteins had matched ids"<<endl;
	  if( toRemove.empty() == false )
	    cout<<"Example of an ID to remove: "<<*toRemove.begin()<<endl;

	  // Remove the matched ids from the training data
	  vGood = filter_ispace( isTrain, toRemove );
	  cout<<"After removing matched ids, training data has ";
	  cout<<isTrain.size()<<" samples; "<<nFeats(isTrain)<<" features"<<endl;
	}
    }

  if( fnAnnots.empty() ) return 0;

  // Loading the Ontology file
  GO::CGOContainer goGraph( fnOBO );

  // Load the GO annotation file
  cout<<"Loading the GO annotations... "; cout.flush();
  GO::CGOACollection goa( fnAnnots );
  cout<<"loaded "<<goa.size()<<" annotations"<<endl;

  // Handle PPI data
  if( fnPPI.empty() == false && toRemove.empty() == false )
    process_ppi( fnPPI, toRemove, sOPrefix, fnIDs );

  // Generate sparse dataset over the output space
  cout<<"Generating output space datasets... "<<endl;
  CDataSet<CSparseSample> dsAnnotMF;
  CDataSet<CSparseSample> dsAnnotBP;
  makeSparseDataset( goa, vGood, dsAnnotMF, goGraph, GO::GO_MF );
  makeSparseDataset( goa, vGood, dsAnnotBP, goGraph, GO::GO_BP );
  cout<<"  MF Annotation data has "<<dsAnnotMF.size()<<" samples; "<<nFeats( dsAnnotMF )<<" features"<<endl;
  cout<<"  BP Annotation data has "<<dsAnnotBP.size()<<" samples; "<<nFeats( dsAnnotBP )<<" features"<<endl;

  // Isolate well-represented samples
  cropSamples( 2, dsAnnotMF );
  cropSamples( 2, dsAnnotBP );
  cout<<"  MF has "<<dsAnnotMF.size()<<" samples; "<<nFeats( dsAnnotMF )<<" features well-represented"<<endl;
  cout<<"  BP has "<<dsAnnotBP.size()<<" samples; "<<nFeats( dsAnnotBP )<<" features well-represented"<<endl;

  // Downsample the training data to the same set of proteins
  vector< string > vFinalMF = dsAnnotMF.getSampleIDs();
  vector< string > vFinalBP = dsAnnotBP.getSampleIDs();

  // Open the files for output
  boost::iostreams::filtering_ostream ofsTrMF;
  boost::iostreams::filtering_ostream ofsTrBP;
  boost::iostreams::filtering_ostream ofsAnnotMF;
  boost::iostreams::filtering_ostream ofsAnnotBP;
  openWriteFile( sOPrefix + "mf_train.sdat.gz", ofsTrMF );
  openWriteFile( sOPrefix + "bp_train.sdat.gz", ofsTrBP );
  openWriteFile( sOPrefix + "mf_annot.sdat.gz", ofsAnnotMF );
  openWriteFile( sOPrefix + "bp_annot.sdat.gz", ofsAnnotBP );

  // Save the results
  isTrain.displaySamples( vFinalMF, ofsTrMF );
  isTrain.displaySamples( vFinalBP, ofsTrBP );
  dsAnnotMF.display( ofsAnnotMF );
  dsAnnotBP.display( ofsAnnotBP );

  return 0;
}
