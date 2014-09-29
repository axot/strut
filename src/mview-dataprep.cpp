// mview-dataprep.cpp - prepares the data for mview experiments
//
// by Artem Sokolov

#include <iostream>
#include <set>
#include <stdarg.h>
#include <fstream>

#include "types.h"
#include "blastout.h"
#include "go-annotation.h"
#include "parsers.h"

using std::cout;
using std::endl;
using std::set;
using std::ofstream;

const GO::ONTOLOGY_INDEX g_myfilter = GO::GO_MF;

typedef CDataSet< CSparseSample > CSparseDataSet;

const string pfxOutput( "mview/yeast/" );

// Finds the set of common IDs
vector< string > commonIDs( const vector< vector< string > >& ids )
{
  unsigned int n = ids.size();

  // Count the occurrences
  sumap_t m;
  for( unsigned int i = 0; i < n; ++i )
    for( vector< string >::const_iterator iter = ids[i].begin();
	 iter != ids[i].end(); ++iter )
      m[ *iter ] = m[ *iter ] + 1;

  // Fetch the intersection
  vector< string > res;
  for( sumap_t::const_iterator iter = m.begin(); iter != m.end(); ++iter )
    if( iter->second == n ) res.push_back( iter->first );

  // Sanity check
  for( unsigned int i = 0; i < n; ++i )
    for( vector< string >::const_iterator iter = res.begin();
	 iter != res.end(); ++iter )
      {
	if( std::find( ids[i].begin(), ids[i].end(), *iter ) == ids[i].end() )
	  throw std::logic_error( "Invalid commonID computation" );
      }

  return res;
}

// Limits a BLAST-hits dataset to examples that contain significant hits in the target organism only
void trimBlastDataSet( CDataSet<CSparseSample>& ds,
		       const vector< string >& targetIDs )
{
  // Traverse the samples
  vector< unsigned int > toKeep;
  for( unsigned int i = 0; i < ds.size(); i++ )
    {
      // Retrieve the sample
      shared_ptr< CSparseSample const > ps = ds.getSample( i );
      
      // Determine if it has any hits in the target organism
      bool bHasHits = false;
      for( unsigned int j = 0; j < targetIDs.size(); j++ )
	{
	  if( ps->getValue( targetIDs[j] ) != 0.0 )
	    bHasHits = true;
	}

      // Keep the sample if it has significant hits
      if( bHasHits == true ) toKeep.push_back( i );
    }

  ds.subsample( toKeep );
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

// Data binder
struct SCrossSpeciesData
{
  SCrossSpeciesData( const GO::CBLASTOutput& blastout )
    : pfmBlast( new CFeatMap ),
      pfmLComplx( new CFeatMap ),
      pfmLocaliz( new CFeatMap ),
      pfmTermini( new CFeatMap ),
      pfmTransmem( new CFeatMap ),
      pfmAnnots( new CFeatMap )
  {
    cout<<"Generating a BLAST dataset... "; cout.flush();
    makeSparseDataset( blastout, dsBlast, 1e-10, 50.0 );
    cout<<dsBlast.size()<<" samples generated"<<endl;
  }

  shared_ptr< CFeatMap > pfmBlast;
  shared_ptr< CFeatMap > pfmLComplx;
  shared_ptr< CFeatMap > pfmLocaliz;
  shared_ptr< CFeatMap > pfmTermini;
  shared_ptr< CFeatMap > pfmTransmem;
  shared_ptr< CFeatMap > pfmAnnots;

  CSparseDataSet dsBlast;
  CSparseDataSet dsLComplx;
  CSparseDataSet dsLocaliz;
  CSparseDataSet dsTermini;
  CSparseDataSet dsTransmem;
  CSparseDataSet dsAnnots;

  void loadI( const string& pfx )
  {
    cout<<"Using prefix: "<<pfx<<endl;

    string fnLComplx = pfx + "-lcomplx.dat.gz";
    cout<<"Loading low-complexity data.... "; cout.flush();
    parseSparseFile( fnLComplx, dsLComplx, pfmLComplx, ',', ' ', ':' );
    cout<<"now at "<<dsLComplx.size()<<" samples, "<<nFeats( dsLComplx )<<" features"<<endl;

    string fnLocaliz = pfx + "-localiz.dat.gz";
    cout<<"Loading localization data.... "; cout.flush();
    parseSparseFile( fnLocaliz, dsLocaliz, pfmLocaliz, ',', ' ', ':' );
    cout<<"now at "<<dsLocaliz.size()<<" samples, "<<nFeats( dsLocaliz )<<" features"<<endl;

    string fnTermini = pfx + "-termini.dat.gz";
    cout<<"Loading termini data.... "; cout.flush();
    parseSparseFile( fnTermini, dsTermini, pfmTermini, ',', ' ', ':' );
    cout<<"now at "<<dsTermini.size()<<" samples, "<<nFeats( dsTermini )<<" features"<<endl;

    string fnTransmem = pfx + "-transmem.dat.gz";
    cout<<"Loading transmembrane data.... "; cout.flush();
    parseSparseFile( fnTransmem, dsTransmem, pfmTransmem, ',', ' ', ':' );
    cout<<"now at "<<dsTransmem.size()<<" samples, "<<nFeats( dsTransmem )<<" features"<<endl;
  }

  void loadO( const string& pfx, const GO::CGOContainer& goGraph )
  {
    string fn = pfx + ".annots.gz";
    cout<<"Loading "<<fn<<"... "; cout.flush();
    GO::CGOACollection goa( fn.c_str() );
    cout<<goa.size()<<" annotations loaded"<<endl;

    cout<<"Generating the annotation samples... "; cout.flush();
    vector< string > ids = IDs();
    makeSparseDataset( goa, ids, dsAnnots, goGraph, g_myfilter, pfmAnnots );
    cout<<"now at "<<dsAnnots.size()<<" samples; "<<nFeats( dsAnnots )<<" features"<<endl;
  }

  void load( const vector< string >& pfxs, const GO::CGOContainer goGraph, const vector< string >& vTargetIDs )
  {
    // Load all input space data first
    for( unsigned int i = 0; i < pfxs.size(); ++i )
      loadI( pfxs[i] );

    // Clip the BLAST hits if necessary
    if( vTargetIDs.empty() == false )
      trimBlast( vTargetIDs );

    // Now the output-space data
    for( unsigned int i = 0; i < pfxs.size(); ++i )
      loadO( pfxs[i], goGraph );
  }

  vector< string > IDs()
  {
    vector< vector< string > > vAll;
    vAll.push_back( dsBlast.getSampleIDs() );    
    vAll.push_back( dsLComplx.getSampleIDs() );  
    vAll.push_back( dsLocaliz.getSampleIDs() );  
    vAll.push_back( dsTermini.getSampleIDs() );  
    vAll.push_back( dsTransmem.getSampleIDs() );
    vector< string > vCommon = commonIDs( vAll );
    return vCommon;
  }

  void trimBlast( const vector< string >& vTargetIDs )
  {
    cout<<"Trimming the BLAST data... "; cout.flush();
    trimBlastDataSet( dsBlast, vTargetIDs );
    cout<<" now at "<<dsBlast.size()<<" samples, "<<nFeats( dsBlast )<<" features"<<endl;
  }

  void subsample( vector< string >& ids )
  {
    dsBlast.subsample( ids );
    dsLComplx.subsample( ids );
    dsLocaliz.subsample( ids );
    dsTermini.subsample( ids );
    dsTransmem.subsample( ids );
    dsAnnots.subsample( ids );
  }

  void save( const string& pfx, vector< string >& ids )
  {
    cout<<"Saving data to prefix "<<pfx<<endl;

    string fn1 = pfx + "-K1.sdat"; ofstream ofs1( fn1.c_str() );
    dsBlast.displaySamples( ids, ofs1 );
    ofs1.close();

    string fn2 = pfx + "-K2.sdat"; ofstream ofs2( fn2.c_str() );
    dsLComplx.displaySamples( ids, ofs2 );
    ofs2.close();

    string fn3 = pfx + "-K3.sdat"; ofstream ofs3( fn3.c_str() );
    dsLocaliz.displaySamples( ids, ofs3 );
    ofs3.close();

    string fn4 = pfx + "-K4.sdat"; ofstream ofs4( fn4.c_str() );
    dsTermini.displaySamples( ids, ofs4 );
    ofs4.close();

    string fn5 = pfx + "-K5.sdat"; ofstream ofs5( fn5.c_str() );
    dsTransmem.displaySamples( ids, ofs5 );
    ofs5.close();

    string fnA = pfx + "-annots.sdat"; ofstream ofsA( fnA.c_str() );
    dsAnnots.displaySamples( ids, ofsA );
  }

  vector< string > featReprO()
  {
    // Count the feature representation in the output space of the target species
    cout<<"Counting the feature representation in the annotation dataset"<<endl;
    simap_t fCount;  nSamplesWFeat( dsAnnots, pfmAnnots, fCount );

    // Find the features that have 10+ samples associated with them
    const int goodCount = 10; vector< string > fGood;
    for( simap_t::const_iterator iter = fCount.begin(); iter != fCount.end(); ++iter )
      if( iter->second >= goodCount ) fGood.push_back( iter->first );
    cout<<"Out of "<<fCount.size()<<" features, "<<fGood.size()<<" are well-represented"<<endl;
    return fGood;
  }

  void remapO( shared_ptr< CFeatMap > pfmA )
  {
    cout<<"Remapping the output space"<<endl;
    remap( dsAnnots, pfmA );
    pfmAnnots = pfmA;
    if( pfmAnnots->nFeats() != nFeats( dsAnnots ) )
      throw std::logic_error( "Output-space remap failed" );
    cout<<"Output space now has "<<dsAnnots.size()<<" samples, "<<nFeats( dsAnnots )<<" features"<<endl;
  }

  void trimO()
  {
    cout<<"Removing annotation samples with fewer than 2 features"<<endl;
    cropSamples( 2, dsAnnots );
    cout<<"Output space now has "<<dsAnnots.size()<<" samples, "<<nFeats( dsAnnots )<<" features"<<endl;
  }
};

struct SWithinSpeciesData
{
  CSparseDataSet dsPPI;
  CSparseDataSet dsAnnots;

  void load( const string& pfx, const GO::CGOContainer& goGraph, shared_ptr< CFeatMap > pfmGood )
  {
    string fnPPI = pfx + "-ppi.sdat.gz";
    cout<<"Loading PPI data.... "; cout.flush();
    parseSparseFile( fnPPI, dsPPI, ',', ',', '=' );
    cout<<" parsed "<<dsPPI.size()<<" samples, "<<nFeats( dsPPI )<<" features"<<endl;
    vector< string > ids = dsPPI.getSampleIDs();

    string fn = pfx + ".annots.gz";
    cout<<"Loading "<<fn<<"... "; cout.flush();
    GO::CGOACollection goa( fn.c_str() );
    cout<<goa.size()<<" annotations loaded"<<endl;

    cout<<"Generating the annotation samples... "; cout.flush();
    makeSparseDataset( goa, ids, dsAnnots, goGraph, g_myfilter );
    cout<<"now at "<<dsAnnots.size()<<" samples; "<<nFeats( dsAnnots )<<" features"<<endl;

    cout<<"Remapping and trimming the output space dataset"<<endl;
    remap( dsAnnots, pfmGood );
    cropSamples( 2, dsAnnots );
    cout<<"Output space is now at "<<dsAnnots.size()<<" samples, "<<nFeats( dsAnnots )<<" features"<<endl;
  }

  void save( const string& pfx, const vector< string >& ids )
  {
    string fnV2I = pfx + "-ppi.sdat"; ofstream ofsV2I( fnV2I.c_str() );
    string fnV2O = pfx + "-annot.sdat"; ofstream ofsV2O( fnV2O.c_str() );
    dsPPI.displaySamples( ids, ofsV2I );
    dsAnnots.displaySamples( ids, ofsV2O );
  }

  void subsample( const vector< string >& ids )
  {
    dsPPI.subsample( ids );
    dsAnnots.subsample( ids );
  }
};

pair< SCrossSpeciesData, SWithinSpeciesData > loadData( const GO::CBLASTOutput& blastout, const GO::CGOContainer& goGraph )
{
  const string pfxTarget = "/s/chopin/c/proj/protfun/data/organisms/s_cerevisiae/2010-02-01/yeast";
  const string pfx1 =      "/s/chopin/c/proj/protfun/data/organisms/d_melanogaster/2009-06-04/fly";
  const string pfx2 =      "/s/chopin/c/proj/protfun/data/organisms/s_pombe/2009-06-04/pombe";

  // Load the target data
  cout<<"---------- Target Species ----------"<<endl;
  SCrossSpeciesData dataTarget( blastout );
  vector< string > pfxsTarget; pfxsTarget.push_back( pfxTarget );
  dataTarget.load( pfxsTarget, goGraph, vector<string>() );

  vector< string > fGood = dataTarget.featReprO();
  shared_ptr< CFeatMap > pfmGood( new CFeatMap( fGood ) );
  dataTarget.remapO( pfmGood );
  dataTarget.trimO();

  vector< string > vTargetIDs = dataTarget.IDs();
  cout<<"There are "<<vTargetIDs.size()<<" target ids"<<endl;

  // Load other species
  cout<<"---------- External Species ----------"<<endl;
  SCrossSpeciesData dataExt( blastout );
  vector< string > pfxsExt; pfxsExt.push_back( pfx1 ); pfxsExt.push_back( pfx2 );
  dataExt.load( pfxsExt, goGraph, vTargetIDs );

  dataExt.remapO( pfmGood );
  dataExt.trimO();

  // Save the external data
  string pfxExtOutput = pfxOutput + "external";
  vector< string > vExtIDs = dataExt.dsAnnots.getSampleIDs();
  dataExt.save( pfxExtOutput, vExtIDs );

  // Load the PPI data
  cout<<"---------- Species-Specific data ----------"<<endl;
  SWithinSpeciesData dataV2;
  dataV2.load( pfxTarget, goGraph, pfmGood );

  // Find the set of common and unique IDs
  vector< string > v1 = dataTarget.dsAnnots.getSampleIDs();
  vector< string > v2 = dataV2.dsAnnots.getSampleIDs();
  std::sort( v1.begin(), v1.end() ); std::sort( v2.begin(), v2.end() );
  vector< string > v12, v1e, v2e;
  std::set_intersection( v1.begin(), v1.end(), v2.begin(), v2.end(),
			 std::inserter( v12, v12.begin() ) );
  std::set_difference( v1.begin(), v1.end(), v2.begin(), v2.end(),
		       std::inserter( v1e, v1e.begin() ) );
  std::set_difference( v2.begin(), v2.end(), v1.begin(), v1.end(),
		       std::inserter( v2e, v2e.begin() ) );

  cout<<"View 1: "<<v1.size()<<endl;
  cout<<"View 2: "<<v2.size()<<endl;
  cout<<"View 1 and 2: "<<v12.size()<<endl;
  cout<<"View 1 exclusively: "<<v1e.size()<<endl;
  cout<<"View 2 exclusively: "<<v2e.size()<<endl;

  if( v1e.size() > 50 ) dataTarget.save( pfxOutput + "v1exclus", v1e );
  if( v2e.size() > 50 ) dataV2.save( pfxOutput + "v2exclus", v2e );

  // Subsample and return the remainder
  std::random_shuffle( v12.begin(), v12.end() );
  dataTarget.subsample( v12 );
  dataV2.subsample( v12 );

  return pair< SCrossSpeciesData, SWithinSpeciesData >( dataTarget, dataV2 );
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

  const string fnBlast = "/s/chopin/c/proj/protfun/data/BLAST/threesp/threesp.blast.gz";
  cout<<"Loading the BLAST hits data... "; cout.flush();
  GO::CBLASTOutput blastout( fnBlast.c_str() );
  cout<<blastout.size()<<" entries loaded"<<endl;

  pair< SCrossSpeciesData, SWithinSpeciesData > data = loadData( blastout, goGraph );

  cout<<"----- Fold Split -----"<<endl;

  cout<<"View 1: "<<endl;
  cout<<"K1: "<<data.first.dsBlast.size()<<" samples, "<<nFeats( data.first.dsBlast )<<" features"<<endl;
  cout<<"K2: "<<data.first.dsLComplx.size()<<" samples, "<<nFeats( data.first.dsLComplx )<<" features"<<endl;
  cout<<"K3: "<<data.first.dsLocaliz.size()<<" samples, "<<nFeats( data.first.dsLocaliz )<<" features"<<endl;
  cout<<"K4: "<<data.first.dsTermini.size()<<" samples, "<<nFeats( data.first.dsTermini )<<" features"<<endl;
  cout<<"K5: "<<data.first.dsTransmem.size()<<" samples, "<<nFeats( data.first.dsTransmem )<<" features"<<endl;
  cout<<"An: "<<data.first.dsAnnots.size()<<" samples, "<<nFeats( data.first.dsAnnots )<<" features"<<endl;
  cout<<endl<<"View 2: "<<endl;
  cout<<"PPI: "<<data.second.dsPPI.size()<<" samples, "<<nFeats( data.second.dsPPI )<<" features"<<endl;
  cout<<"Ann: "<<data.second.dsAnnots.size()<<" samples, "<<nFeats( data.second.dsAnnots )<<" features"<<endl;

  vector< string > ids = data.first.dsAnnots.getSampleIDs();

  // Find connected components in the thresholded e-value graph
  const double pi_thresh = 40;
  vector< vector< string > > conComp;
  for( unsigned int i = 0; i < ids.size(); ++i )
    {
      const string id = ids[i];

      // Display progress
      if( i % 100 == 0 ) {cout<<"."; cout.flush();}

      // Traverse the existing connected components and find nearby ones
      vector< unsigned int > vi;
      for( unsigned int j = 0; j < conComp.size(); j++ )
	for( unsigned int k = 0; k < conComp[j].size(); k++ )
	  {
	    //	    if( data.blastout->proximityEVal( id, conComp[j][k], e_thresh ) == true )
	    if( blastout.proximityPIden( id, conComp[j][k], pi_thresh ) == true )
	      {
		vi.push_back( j );
		break;		// for( k = 0; ...
	      }
	  }

      // Add a new component
      if( vi.size() == 0 )
	{
	  vi.push_back( conComp.size() );
	  vector< string > v;
	  conComp.push_back( v );
	}

      unsigned int j0 = vi[0];

      // Merge components as needed
      if( vi.size() > 1 )
	{
	  for( unsigned int jj = 1; jj < vi.size(); jj++ )
	    {
	      unsigned int j = vi[jj];

	      // Move the contents of the j^th component to j0^th
	      for( unsigned int k = 0; k < conComp[j].size(); k++ )
		conComp[j0].push_back( conComp[j][k] );

	      // Clear the j^th component
	      conComp[j].clear();
	    }
	}

      // Append the current id
      conComp[ j0 ].push_back( id );
    }
  cout<<endl;

  cout<<"Sizes of conComp: ";
  for( unsigned int i = 0; i < conComp.size(); i++ )
    cout<<conComp[i].size()<<" ";
  cout<<endl;

  // Additional checks
  
  // Every protein must be within thresh of something else in a component
  cout<<"Checking within-component constraints"<<endl;
  for( unsigned int i = 0; i < conComp.size(); i++ )
    {
      if( conComp[i].size() < 2 ) continue;
      for( unsigned int j = 0; j < conComp[i].size(); j++ )
      {
	bool bHasNeigh = false;
	for( unsigned int k = 0; k < conComp[i].size(); k++ )
	  {
	    if( j == k ) continue;
	    bHasNeigh |= blastout.proximityPIden( conComp[i][j], conComp[i][k], pi_thresh );
	    //	    bHasNeigh |= data.blastout->proximityEVal( conComp[i][j], conComp[i][k], e_thresh );
	  }

	if( bHasNeigh == false )
	  {
	    cout<<"Component "<<i<<endl;
	    throw std::logic_error( "A component is not fully connected" );
	  }
      }
    }

  // There must be no proteins within thresh of each other across components
  cout<<"Checking cross-component constraints"<<endl;
  for( unsigned int i = 0; i < conComp.size(); i++ )
    for( unsigned int j = i+1; j < conComp.size(); j++ )
      for( unsigned int ki = 0; ki < conComp[i].size(); ki++ )
	for( unsigned int kj = 0; kj < conComp[j].size(); kj++ )
	  {
	    if(blastout.proximityPIden(conComp[i][ki], conComp[j][kj], pi_thresh))
	      //	    if(data.blastout->proximityEVal(conComp[i][ki], conComp[j][kj], e_thresh))
	      {
		cout<<"Components "<<i<<" and "<<j<<endl;
		throw std::logic_error( "Components are connected" );
	      }
	  }
  
  // Count the number of non-zero components and find the largest one
  unsigned int nNZ = 0;
  unsigned int iLargest = 0;
  unsigned int nLargest = 0;
  for( unsigned int i = 0; i < conComp.size(); i++ )
    {
      if( conComp[i].size() > 0 ) nNZ++;
      if( conComp[i].size() > nLargest )
	{
	  nLargest = conComp[i].size();
	  iLargest = i;
	}
    }

  cout<<"Number of connected components: "<<nNZ<<endl;
  cout<<"The size of the largest component: "<<nLargest<<endl;

  // Compose the fold split
  const unsigned int nFolds = 5;
  vector< vector< string > > folds( nFolds );
  
  while( nNZ > 0 )
    {
      // Find the smallest fold
      unsigned int iFold = 0;
      unsigned int iSize = folds[0].size();
      for( unsigned int i = 0; i < nFolds; i++ )
	if( folds[i].size() < iSize )
	  {
	    iSize = folds[i].size();
	    iFold = i;
	  }

      // Move the data from the largest component to the smallest fold
      for( unsigned int i = 0; i < conComp[iLargest].size(); i++ )
	folds[iFold].push_back( conComp[iLargest][i] );
      conComp[iLargest].clear();

      // Recompute the number of non-zero components and the largest component
      nNZ = 0; iLargest = 0; nLargest = 0;
      for( unsigned int i = 0; i < conComp.size(); i++ )
	{
	  if( conComp[i].size() > 0 ) nNZ++;
	  if( conComp[i].size() > nLargest )
	    {
	      nLargest = conComp[i].size();
	      iLargest = i;
	    }
	}
    }

  cout<<"Optimal fold split: ";
  for( unsigned int i = 0; i < nFolds; i++ )
    cout<<folds[i].size()<<" ";
  cout<<endl;

  // Additional check
  cout<<"Checking cross-fold constraints"<<endl;
  for( unsigned int i = 0; i < nFolds; i++ )
    for( unsigned int j = i+1; j < nFolds; j++ )
      for( unsigned int ki = 0; ki < folds[i].size(); ki++ )
	for( unsigned int kj = 0; kj < folds[j].size(); kj++ )
	  {
	    if( folds[i][ki] == folds[j][kj] )
	      {
		cout<<"Folds "<<i<<" and "<<j<<" : "<<folds[i][ki]<<endl;
		throw std::logic_error( "Protein present in multiple folds" );
	      }

	    if(blastout.proximityPIden(folds[i][ki], folds[j][kj], pi_thresh))
	      {
		cout<<"Folds "<<i<<" and "<<j<<endl;
		throw std::logic_error( "Folds violate the proximity constraints" );
	      }
	  }

  // Save the fold data
  for( unsigned int i = 0; i < nFolds; i++ )
    {
      // Setup the filenames
      string sIndex = boost::lexical_cast< string >( i );
      string fn1 = pfxOutput + "v1-fold" + sIndex;
      string fn2 = pfxOutput + "v2-fold" + sIndex;

      // Save the data
      data.first.save( fn1, folds[i] );
      data.second.save( fn2, folds[i] );

      // Compose the filenames
      // string sIndex = boost::lexical_cast< string >( i );
      // string fnB = fnPrefix + "blast_sc_" + sIndex + ".sdat";
      // string fnP = fnPrefix + "ppi_sc_" + sIndex + ".sdat";
      // //      string fnA = fnPrefix + "annot_sc_" + sIndex + ".sdat";

      // // Open the streams
      // std::ofstream ofB( fnB.c_str() );
      // std::ofstream ofP( fnP.c_str() );
      // //      std::ofstream ofA( fnA.c_str() );

      // // Save the data
      // data.pdsBlast->displaySamples( folds[i], ofB );
      // data.pdsPPI->displaySamples( folds[i], ofP );
      // //      data.pdsAnnotSC->displaySamples( folds[i], ofA );

      // // Close the streams
      // ofB.close(); ofP.close(); //ofA.close();
    }
  
  return 0;
}

