// gostruct.cpp - prediction of GO terms using structured-outputs framework
//
// by Artem Sokolov

#include "misc.h"
#include "eval.h"
#include "parsers.h"

#include <iostream>
#include <string>
#include <fstream>

#include <map>

#include <memory>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::multimap;
using std::pair;
using std::shared_ptr;

void findSRate( const CDataSet<CSparseSample>& ds, const CDataSet<CSparseSample>& dsTruth,
		map<string, double>& thresholds )
{
  // Retrieve the feature maps
  shared_ptr< CFeatMap const > pfm = ds.getSample( 0 )->getFeatMap();
  shared_ptr< CFeatMap const > pfmTruth = ds.getSample( 0 )->getFeatMap();

  for( unsigned int iFeat = 0; iFeat < nFeats( ds ); iFeat++ )
    {
      string s = pfm->i2f( iFeat );
      int jFeat = pfmTruth->f2i( s );
      if( jFeat < 0 ) throw std::logic_error( "Mismatched datasets" );

      multimap< double, bool > mymap;

      // Pull out the labels corresponding to this feature
      for( unsigned int i = 0; i < ds.size(); i++ )
	{
	  string ss = ds.i2s( i );
	  int j = dsTruth.s2i( ss );
	  if( j < 0 ) throw std::logic_error( "Mismatched datasets" );

	  pair< double, bool > entry;
	  entry.first = ds.getSample(i)->getValue( iFeat );
	  entry.second = dsTruth.getSample(j)->getValue( jFeat ) > 0.0;

	  mymap.insert( entry );
	}

      // Find the corresponding threshold
      pair< string, double > p;
      p.first = s;
      p.second = bestSRate( mymap );
      thresholds.insert( p );
    }
}

// Computes the average number of non-zero features in a dataset
// Computes the min and the max as well
double aveNZFeats( const CDataSet<CSparseSample>& ds, unsigned int& nMin, unsigned int& nMax )
{
  double res = 0.0; nMin = std::numeric_limits<unsigned int>::max(); nMax = 0;
  for( unsigned int i = 0; i < ds.size(); i++ )
    {
      unsigned int k = ds.getSample(i)->L0();
      res += k;
      if( k < nMin ) nMin = k;
      if( k > nMax ) nMax = k;
    }
  return res / static_cast<double>(ds.size());
}

// Thresholds a dataset at the top-scoring k samples for each node
// Returns a copy
shared_ptr< CDataSet<CSparseSample> > threshDataSet( const CDataSet<CSparseSample>& ds, map< string, unsigned int > kmap )
{
  // Create a copy
  shared_ptr< CDataSet<CSparseSample> > pRes( new CDataSet<CSparseSample>( ds ) );
  threshold( *pRes, kmap );

  return pRes;
}

// Thresholds a dataset using a double value for each node
// Returns a copy
shared_ptr< CDataSet<CSparseSample> > threshDataSet( const CDataSet<CSparseSample>& ds, map< string, double > prof )
{
  // Create a copy
  shared_ptr< CDataSet<CSparseSample> > pRes( new CDataSet<CSparseSample>( ds ) );
  threshold( *pRes, prof );

  return pRes;
}

// Takes truth and predicted datasets and computes kernel loss, precision, recall and ave. #nodes
void computeAll( const CDataSet<CSparseSample>& dsTruth, const CDataSet<CSparseSample> dsPred,
		 double& loss, double& precision, double& recall, double& nNodes )
{
  loss = dsPred.loss( dsTruth );
  pair<double, double> pnr = computePnR( dsPred, dsTruth );
  precision = pnr.first;
  recall = pnr.second;

  unsigned int min, max;
  nNodes = aveNZFeats( dsPred, min, max );
}

// Computes and displays loss, precision, recall and the number of nodes
void doWork( const CDataSet<CSparseSample>& ds1, const CDataSet<CSparseSample>& ds2, string name )
{
  double loss, precision, recall, nNodes;
  computeAll( ds1, ds2, loss, precision, recall, nNodes );
  cout<<name<<": loss = "<<loss<<"  prc = "<<precision<<"  rec = "<<recall<<"  #n = "<<nNodes<<endl;
}

// Entry point of the program
int main( int argc, char* argv[] )
{
  if( argc < 2 )
    {
      cout<<"Usage: "<<argv[0]<<" <truth matrix> <optional pred matrix> <optional output thresh matrix>"<<std::endl;
      return -1;
    }

  // Create a kernel and a loss
  typedef CKernelLoss<CSparseSample> CSparseKLoss;
  CDataSet<CSparseSample>::binop_t fker = CSparseKernel(true);
  CDataSet<CSparseSample>::binop_t floss = CSparseKLoss(fker);

  // Create the datasets
  CDataSet<CSparseSample> dsTruth( fker, floss );
  CDataSet<CSparseSample> dsPred( fker, floss );

  // Load the first dataset
  parseTabDelFile( argv[1], dsTruth );

  // Perform truth matrix analysis
  if( argc == 2 )
    {
      for( unsigned int i = 0; i < dsTruth.size(); i++ )
	cout<<dsTruth.getSample( i )->L0()<<endl;
      return 0;
    }

  cout<<"Truth: "<<argv[1]<<" has "<<dsTruth.size()<<" samples; "<<nFeats( dsTruth )<<" features"<<std::endl;

  // Load the second dataset
  string fnPred( argv[2] );
  if( fnPred.find( ".tdel" ) != fnPred.npos )
    parseTabDelFile( argv[2], dsPred );
  else if( fnPred.find( ".sdat" ) != fnPred.npos )
    parseSparseFile( argv[2], dsPred, ',', ',', '=' );
  else
    throw std::runtime_error( "Unable to determine the format" );

  cout<<"Predictions: "<<argv[2]<<" has "<<dsPred.size()<<" samples; "<<nFeats( dsPred )<<" features"<<std::endl;

  // Determine the intersection of samples/features between the two datasets
  vector< string > vCommonSmpls = commonSampleIDs<CSparseSample>( dsTruth, dsPred );
  vector< string > vCommonFeats = commonFeatIDs( dsTruth, dsPred );
  cout<<"The intersection has "<<vCommonSmpls.size()<<" samples; "<<
    vCommonFeats.size()<<" features"<<endl;

  // Create the feature map for the common features
  shared_ptr< CFeatMap > pfm( new CFeatMap( vCommonFeats ) );

  // Trim both datasets to the intersection
  dsTruth.subsample( vCommonSmpls ); remap( dsTruth, pfm );
  dsPred.subsample( vCommonSmpls ); remap( dsPred, pfm );
  cout<<"After trimming truth has "<<dsTruth.size()<<" samples; "<<nFeats( dsTruth )<<" features"<<std::endl;
  cout<<"After trimming predictions have "<<dsPred.size()<<" samples; "<<nFeats( dsPred )<<" features"<<std::endl;

  // Compute the min, max and average statistics
  unsigned int min1, max1;
  unsigned int min2, max2;
  double ave1 = aveNZFeats(dsTruth, min1, max1);
  double ave2 = aveNZFeats(dsPred, min2, max2);

  // Display some basic statistic about the matrices
  cout<<"Truth - number of annotations: min = "<<min1<<"  max = "<<max1<<"  ave = "<<ave1<<endl;
  cout<<"Preds - number of annotations: min = "<<min2<<"  max = "<<max2<<"  ave = "<<ave2<<endl;

  // Retrieve the range of predictions
  std::pair<double, double> pRange = getRange( dsPred );
  cout<<"Predictions are in ( "<<pRange.first<<", "<<pRange.second<<" )"<<std::endl;

  // Generate and save the ROC curves
  string fnROC = string( argv[2] ) + ".roc";
  cout<<"Saving the ROC curves to "<<fnROC<<endl;
  ofstream ofsROC( fnROC.c_str() );
  double auroc = 0.0;
  for( unsigned int fi = 0; fi < pfm->nFeats(); ++fi )
    {
      string f = pfm->i2f( fi );

      // Compose the (score, label) pairs
      vector< pair< double, unsigned int > > rocData;
      for( unsigned int i = 0; i < vCommonSmpls.size(); ++i )
	{
	  if( dsTruth.i2s( i ) != dsPred.i2s( i ) ) throw std::logic_error( "Inconsisted data" );

	  double score = dsPred.getSample( i )->getValue( f );
	  unsigned int label = dsTruth.getSample( i )->getValue( f ) != 0.0;
	  pair< double, unsigned int > p( score, label );
	  rocData.push_back( p );
	}

      // Compute the curve
      vector< pair< double, double > > roc = ROC( rocData );
      auroc += AUROC( roc );
      
      // Save the results to a file
      ofsROC<<f<<" FP"; for( unsigned int i = 0; i < roc.size(); ++i ) ofsROC<<" "<<roc[i].first; ofsROC<<endl;
      ofsROC<<f<<" TP"; for( unsigned int i = 0; i < roc.size(); ++i ) ofsROC<<" "<<roc[i].second; ofsROC<<endl;
    }
  ofsROC.close();
  auroc /= static_cast<double>( pfm->nFeats() );
  cout<<"Average area under the ROC: "<<auroc<<endl;

  // No thresholding
  doWork( dsTruth, dsPred, "Raw" );

  if( argc == 4 )
    {
      // Load the profile dataset
      shared_ptr< CDataSet<CSparseSample> > pProf( new CDataSet<CSparseSample>(fker, floss) );
      parseTabDelFile( argv[3], *pProf );
      cout<<"Profile: "<<argv[3]<<" has "<<pProf->size()<<" samples; "<<nFeats( pProf )<<" features"<<endl;

      // Trim the profile
      pProf->subsample( vCommonSmpls ); remap( *pProf, pfm );
      if( pProf->size() != dsPred.size() || nFeats( pProf ) != nFeats( dsPred ) )
	throw std::logic_error( "Profile does not match the set of predictions" );
      cout<<"After trimming Profile has "<<pProf->size()<<" samples; "<<nFeats( pProf )<<" features"<<endl;

      // Extract the profile
      map< string, unsigned int > prof;
      for( unsigned int iFeat = 0; iFeat < nFeats( pProf ); iFeat++ )
	{
	  string f = pfm->i2f( iFeat );
	  prof[f] = nSamplesWFeat( *pProf, f );
	}
      
      // Threshold based on the profile
      shared_ptr< CDataSet<CSparseSample> > pThresh = threshDataSet( dsPred, prof );
      doWork( dsTruth, *pThresh, "prof" );
    }
  else
    {
      // Compute the success rate thresholds
      map< string, double > thresholds;
      findSRate( dsPred, dsTruth, thresholds );

      // Threshold the dataset
      shared_ptr< CDataSet<CSparseSample> > pThresh = threshDataSet( dsPred, thresholds );
      doWork( dsTruth, *pThresh, "BSR" );
    }

  

  return 0;

}
