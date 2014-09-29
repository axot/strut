// mview-dataprep-mouse.cpp - prepares the mouse data for mview experiments
//
// by Artem Sokolov

#include <iostream>

#include "types.h"
#include "parsers.h"
#include "blastout.h"
#include "go-annotation.h"

using std::cout;
using std::endl;

const GO::ONTOLOGY_INDEX g_myfilter = GO::GO_MF;

const string blast_filename( "/s/chopin/c/proj/protfun/data/BLAST/human_and_mouse/human_and_mouse.blast.gz" );
const string ppi_filename( "/s/chopin/c/proj/protfun/data/organisms/m_musculus/mouse_ppi.sdat" );
const string humanAnnot_filename( "/s/chopin/c/proj/protfun/data/organisms/h_sapiens/gene_association.goa_human" );
const string mouseAnnot_filename( "/s/chopin/c/proj/protfun/data/organisms/m_musculus/gene_association.mgi" );

const string fnPrefix( "mview/m_musculus/mf_" );

// Set intersection between two sorted vectors with an intuitive interface
template< typename _T >
vector<_T> set_intersection( const vector<_T>& u, const vector<_T>& v )
{
  vector<_T> res;
  std::set_intersection( u.begin(), u.end(), v.begin(), v.end(),
			 std::inserter( res, res.begin() ) );
  return res;
}

// Set difference between two sorted vectors with an intuitive interface
template< typename _T >
vector<_T> set_difference( const vector<_T>& u, const vector<_T>& v )
{
  vector<_T> res;
  std::set_difference( u.begin(), u.end(), v.begin(), v.end(),
		       std::inserter( res, res.begin() ) );
  return res;
}

// Limits a BLAST-hits dataset to examples that contain significant hits in the target organism only
void trimBlastDataSet( shared_ptr< CDataSet<CSparseSample> > pds,
		       const vector< string >& targetIDs )
{
  // Traverse the samples
  vector< unsigned int > toKeep;
  for( unsigned int i = 0; i < pds->size(); i++ )
    {
      // Retrieve the sample
      shared_ptr< CSparseSample const > ps = pds->getSample( i );
      
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

  pds->subsample( toKeep );
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

  cout<<"Generating a BLAST dataset... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsBlast( new CDataSet<CSparseSample> );
  makeSparseDataset( blastout, *pdsBlast, 1e-10, 50.0 );
  cout<<pdsBlast->size()<<" samples generated"<<endl;

  cout<<"Loading PPI data... "; cout.flush();
  shared_ptr< CDataSet<CSparseSample> > pdsPPI( new CDataSet<CSparseSample> );
  parseSparseFile( ppi_filename, *pdsPPI, ',', ',', '=' );
  cout<<pdsPPI->size()<<" samples loaded"<<endl;

  cout<<"Loading human annotations... "; cout.flush();
  GO::CGOACollection goaHuman( humanAnnot_filename.c_str() );
  cout<<goaHuman.size()<<" annotations loaded"<<endl;

  cout<<"Loading mouse annotations... "; cout.flush();
  GO::CGOACollection goaMouse( mouseAnnot_filename.c_str() );
  cout<<goaMouse.size()<<" annotations loaded"<<endl;

  // Set up the IDs and feature map for datasets
  vector< string > sidsBlast = pdsBlast->getSampleIDs();
  shared_ptr< CFeatMap > pfm( new CFeatMap );
  shared_ptr< CDataSet<CSparseSample> > pdsHumanGOA( new CDataSet<CSparseSample> );
  shared_ptr< CDataSet<CSparseSample> > pdsMouseGOA( new CDataSet<CSparseSample> );

  cout<<"Generating human annotation dataset... "; cout.flush();
  makeSparseDataset( goaHuman, sidsBlast, *(pdsHumanGOA), goGraph, g_myfilter, pfm );
  cout<<pdsHumanGOA->size()<<" samples generated"<<endl;

  cout<<"Generating mouse annotation dataset... "; cout.flush();
  makeSparseDataset( goaMouse, sidsBlast, *(pdsMouseGOA), goGraph, g_myfilter, pfm );
  cout<<pdsMouseGOA->size()<<" samples generated"<<endl;

  cout<<"Trimming the BLAST set to samples with significant hits in the target organism"<<endl;
  cout<<"Before trim: "<<pdsBlast->size()<<" samples; "<<nFeats( pdsBlast )<<" features"<<endl;
  trimBlastDataSet( pdsBlast, pdsMouseGOA->getSampleIDs() );
  cout<<"After trim: "<<pdsBlast->size()<<" samples; "<<nFeats( pdsBlast )<<" features"<<endl;
  
  // Trim the annotations accordingly
  sidsBlast = pdsBlast->getSampleIDs();
  pdsHumanGOA->subsample( sidsBlast, false );
  pdsMouseGOA->subsample( sidsBlast, false );
  cout<<"Human annotation dataset now has "<<pdsHumanGOA->size()<<" samples"<<endl;
  cout<<"Mouse annotation dataset now has "<<pdsMouseGOA->size()<<" samples"<<endl;

  // Count the number of mouse examples associated with each feature
  simap_t featCount;
  for( unsigned int j = 0; j < pfm->nFeats(); j++ )
    {
      unsigned int n = 0;
      for( unsigned int i = 0; i < pdsPPI->size(); i++ )
	if( pdsPPI->getSample( i )->getValue( j ) != 0.0 ) n++;

      string f = pfm->i2f( j );
      featCount[f] = n;
    }

  // Project onto the features that are represented in at least 10 examples
  vector< string > reprFeat;
  for( simap_t::const_iterator iter = featCount.begin(); iter != featCount.end(); iter++ )
    if( iter->second >= 10 ) reprFeat.push_back( iter->first );
  shared_ptr< CFeatMap > pfmRepr( new CFeatMap( reprFeat ) );
  remap( *pdsHumanGOA, pfmRepr );
  remap( *pdsMouseGOA, pfmRepr );
  cout<<"Out of "<<featCount.size()<<" features, "<<pfmRepr->nFeats()<<" are well-represented"<<endl;

  // Remove samples that contain no features or root only
  cout<<"Removing annotation sample with fewer than 2 features"<<endl;
  cropSamples( 2, *pdsHumanGOA );
  cropSamples( 2, *pdsMouseGOA );
  cout<<"Human annotation dataset now has "<<pdsHumanGOA->size()<<" samples"<<endl;
  cout<<"Mouse annotation dataset now has "<<pdsMouseGOA->size()<<" samples"<<endl;

  // Retrieve the remaining set of ids
  vector< string > idsHuman = pdsHumanGOA->getSampleIDs();
  vector< string > idsMouse = pdsMouseGOA->getSampleIDs();

  // Save the human data and mouse annotations
  string fnHumanBlast = fnPrefix + "human_blast.sdat";
  string fnHumanGOA = fnPrefix + "human_annot.sdat";
  string fnMouseGOA = fnPrefix + "mouse_annot.sdat";
  std::ofstream ofHumanBlast( fnHumanBlast.c_str() );
  //  std::ofstream ofHumanGOA( fnHumanGOA.c_str() );
  //  std::ofstream ofMouseGOA( fnMouseGOA.c_str() );
  pdsBlast->CDataSet<CSparseSample>::displaySamples( idsHuman, ofHumanBlast );
  pdsHumanGOA->save( fnHumanGOA );
  pdsMouseGOA->save( fnMouseGOA );
  //  pdsHumanGOA->CDataSet<CSparseSample>::displaySamples( idsHuman, ofHumanGOA );
  //  pdsMouseGOA->CDataSet<CSparseSample>::displaySamples( idsMouse, ofMouseGOA );

  cout<<"Downsampling to the mouse data of interest..."<<endl;
  pdsPPI->subsample( idsMouse, false );
  pdsBlast->subsample( idsMouse );
  if(  pdsBlast->size() != pdsMouseGOA->size() )
    throw std::logic_error( "Mismatched datasets" );
  
  sidsBlast = pdsBlast->getSampleIDs();
  vector< string > sidsPPI = pdsPPI->getSampleIDs();
  vector< string > sidsCommon = set_intersection( sidsBlast, sidsPPI );
  vector< string > sidsBLASTe = set_difference( sidsBlast, sidsCommon );
  vector< string > sidsPPIe = set_difference( sidsPPI, sidsCommon );

  cout<<"# BLAST  ids : "<<sidsBlast.size()<<endl;
  cout<<"# PPI    ids : "<<sidsPPI.size()<<endl;
  cout<<"# common ids : "<<sidsCommon.size()<<endl;
  cout<<"# BLASTe ids : "<<sidsBLASTe.size()<<endl;
  cout<<"# PPIe   ids : "<<sidsPPIe.size()<<endl;

  // Save mouse BLAST samples that are missing matching PPI data
  string fnMouseBlast = fnPrefix + "mouse_blast_noppi.sdat";
  std::ofstream ofMouseBlast( fnMouseBlast.c_str() );
  pdsBlast->CDataSet<CSparseSample>::displaySamples( sidsBLASTe, ofMouseBlast );

  cout<<"----------------------"<<endl;
  cout<<"Performing fold split"<<endl;
  const double pi_thresh = 50;

  // Find connected components in the thresholded percent-identity graph
  vector< vector< string > > conComp;
  for( unsigned int i = 0; i < sidsCommon.size(); i++ )
    {
      const string id = sidsCommon[i];
      
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
      // Compose the filenames
      string sIndex = boost::lexical_cast< string >( i );
      string fnB = fnPrefix + "mouse_blast_" + sIndex + ".sdat";
      string fnP = fnPrefix + "mouse_ppi_" + sIndex + ".sdat";

      // Open the streams
      std::ofstream ofB( fnB.c_str() );
      std::ofstream ofP( fnP.c_str() );

      // Save the data
      pdsBlast->displaySamples( folds[i], ofB );
      pdsPPI->displaySamples( folds[i], ofP );

      // Close the streams
      ofB.close(); ofP.close();
    }
  
  return 0;
}
