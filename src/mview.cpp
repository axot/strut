// mview.cpp - multi-view experiments
//
// by Artem Sokolov

#include "sample.h"
#include "kernel.h"
#include "loss.h"
#include "io-dataset.h"
#include "clsf.h"
#include "types.h"
#include "cosvm.h"
#include "parsers.h"
#include "params.h"
#include "eval.h"

#include <numeric>

#include <iostream>

#include <boost/lexical_cast.hpp>

using std::cout;
using std::endl;
using std::set;

const string pfxV1Pred = "/s/chopin/c/proj/protfun/users/sokolov/mview/pred/yeast/exp1/view1";

// Checks for duplicate sample IDs and throws an exception if any found
template< typename _T >
void checkDuplicates( shared_ptr< CDataSet<_T> const > pds )
{
  vector< string > ids = pds->getSampleIDs();
  set< string > uids;
  for( unsigned int i = 0; i < ids.size(); i++ )
    uids.insert( ids[i] );
  if( ids.size() != uids.size() )
    throw std::logic_error( "Dataset contains duplicates" );
}

// A bundle data strucure
struct SData
{
  typedef CIODataSet<vSparseSample, CSparseSample> CIOSet;
  shared_ptr< CIOSet > pTrainV1;
  shared_ptr< CIOSet > pTrainV2;
  shared_ptr< CIOSet > pTestV1;
  shared_ptr< CIOSet > pTestV2;

  SData( const CDataSet<vSparseSample>::binop_t& fiker1,
	 const CDataSet<vSparseSample>::binop_t& fiker2,
	 shared_ptr< CDataSet<CSparseSample> > pO,
	 const function< double(double, double) >& fioker )
    : pTrainV1( new CIOSet( fiker1, pO, fioker ) ),
      pTrainV2( new CIOSet( fiker2, pO, fioker ) ),
      pTestV1( new CIOSet( fiker1, pO, fioker ) ),
      pTestV2( new CIOSet( fiker2, pO, fioker ) ) {}

  void display()
  {
    cout<<"View1: Train - "<<pTrainV1->sizeI()<<" samples; Test - "<<pTestV1->sizeI()<<" samples; ";
    cout<<nKernels( *pTrainV1->getI() )<<" kernels for "<<nFeats( *pTrainV1->getI() )<<" features total"<<endl;
    cout<<"View2: Train - "<<pTrainV2->sizeI()<<" samples; Test - "<<pTestV2->sizeI()<<" samples; ";
    cout<<nKernels( *pTrainV2->getI() )<<" kernels for "<<nFeats( *pTrainV2->getI() )<<" features total"<<endl;
    cout<<"Output space: "<<pTrainV1->sizeO()<<" samples; "<<nFeats( *pTrainV1->getO() )<<" features"<<endl;
  }
};

// Loads a single kernel space for all folds
// Returns the feature map associated with the space
shared_ptr< CFeatMap > loadFoldSpace( const string& pfx, const string& sfx,
				      unsigned int nFolds,
				      unsigned int testFold,
				      CDataSet< CSparseSample >& dsTr,
				      CDataSet< CSparseSample >& dsTe )
{
  // Setup a common feature map
  shared_ptr< CFeatMap > pfm( new CFeatMap );

  // Load the common data
  for( unsigned int i = 0; i < nFolds; i++ )
    {
      // Compose the filenames
      string sIndex = boost::lexical_cast< string >( i );
      string fn = pfx + "-fold" + sIndex + sfx;
      //      cout<<"Loading "<<fn<<endl;

      // Handle the test fold
      if( i == testFold )
	parseSparseFile( fn, dsTe, pfm, ',', ',', '=' );
      else
	parseSparseFile( fn, dsTr, pfm, ',', ',', '=' );
    }  

  return pfm;
}

// Loads the view 1 prediction data
void loadV1Preds( const string& pfx, const string& sfx, 
		  unsigned int nFolds, unsigned int testFold,
		  CDataSet<CSparseSample>& dsTr,
		  CDataSet<CSparseSample>& dsTe )
{
  // Setup a common feature map
  shared_ptr< CFeatMap > pfm( new CFeatMap );

  for( unsigned int i = 0; i < nFolds; ++i )
    {
      // Compose the filename that contains View 1 predictions
      string fn = pfx + boost::lexical_cast<string>( i ) + sfx;
      cout<<"Loading View 1 predictions from "<<fn<<endl;

      if( i == testFold )
	parseSparseFile( fn, dsTe, pfm, ',', ',', '=' );
      else
	parseSparseFile( fn, dsTr, pfm, ',', ',', '=' );
    }
}

// Performes sample number adjustment in place
void adjustNSamples( unsigned int nAvailable, unsigned int& nLabeled, unsigned int& nUnlabeled )
{
  cout<<"Number of requested labeled training samples: "<<nLabeled<<endl;
  cout<<"Number of requested unlabeled training samples: "<<nUnlabeled<<endl;
  cout<<"Number of available samples: "<<nAvailable<<endl;
  
  // Adjust the number of labeled and unlabeled examples
  if( nLabeled > nAvailable )
    {
      nLabeled = nAvailable;
      nUnlabeled = 0;
    }
  else if( nLabeled + nUnlabeled > nAvailable )
  nUnlabeled = nAvailable - nLabeled;

  cout<<"Number of labeled training examples to be used: "<<nLabeled<<endl;
  cout<<"Number of unlabeled training examples to be used: "<<nUnlabeled<<endl;
}

// Appends data to the IO dataset
template< typename _T >
void addData( CIODataSet<_T, CSparseSample>& iods,
	      const CDataSet<_T>& ids, const CDataSet<CSparseSample>& ods,
	      unsigned int nLabeled, unsigned int nUnlabeled )
{
  unsigned int nTotal = nLabeled + nUnlabeled;
  if( ids.size() != nTotal ) throw std::logic_error( "Bad arguments to addData()" );

  // Add unlabeled training examples
  for( unsigned int i = nLabeled; i < nTotal; ++i )
    {
      string s = ids.i2s( i );
      shared_ptr< const _T > p = ids.getSample( i );
      iods.addInputSample( s, p, std::numeric_limits<unsigned int>::max() );
    }

  // Add labeled training examples
  for( unsigned int i = 0; i < nLabeled; ++i )
    {
      string s = ids.i2s( i );
      int j = ods.s2i( s );
      if( j < 0 ) throw std::logic_error( "Missing annotation for " + s );
      
      shared_ptr< const _T > pi = ids.getSample( i );
      shared_ptr< const CSparseSample > po = ods.getSample( j );
      iods.addSample( s, pi, po );
    }
}

// Loads all the data
SData loadData( const string& pfx, unsigned int nFolds, unsigned int testFold,
		const CStrutParams& params, unsigned int nLabeled, unsigned int nUnlabeled )
{
  const string pfxV1P = pfxV1Pred + "-" + boost::lexical_cast<string>( nLabeled ) + "/";
  const string fnExt = pfx + "external";
  const string fnV1 = pfx + "v1";
  const string fnV2 = pfx + "v2";

  // Composite datasets
  CDataSet< vSparseSample > dsTrV1, dsTeV1, dsExV1;
  CDataSet< vSparseSample > dsTrV2, dsTeV2;
  CDataSet< CSparseSample > dsTrPPI, dsTePPI;
  CDataSet< CSparseSample > dsAnnots;
  CDataSet< CSparseSample > dsDud;

  // Unlabeled examples are only supported by the "co" and "trans" algorithms
  if( params.alg_choice() != "co" && params.alg_choice() != "trans" )
    nUnlabeled = 0;

  // Load the PPI data for view 2
  const string sfxPPI = "-ppi.sdat";
  loadFoldSpace( fnV2, sfxPPI, nFolds, testFold, dsTrPPI, dsTePPI );

  // Fetch the number of available samples and perform the dataset size adjustment
  adjustNSamples( dsTrPPI.size(), nLabeled, nUnlabeled );

  // Downsample the PPI data
  unsigned int nTotal = nLabeled + nUnlabeled;
  vector< unsigned int > v( nTotal );
  for( unsigned int i = 0; i < nTotal; ++i ) v[i] = i;
  dsTrPPI.subsample( v );

  // Load the kernel spaces for view 1
  const unsigned int nKernels = 5;
  for( unsigned int i = 1; i <= nKernels; ++i )
    {
      CDataSet< CSparseSample > dsTrK, dsTeK, dsExK;
      
      // Compose the filename suffix
      const string sfx = "-K" + boost::lexical_cast<string>(i) + ".sdat";

      // Load the data for the current kernel
      shared_ptr< CFeatMap > pfm = loadFoldSpace( fnV1, sfx, nFolds, testFold, dsTrK, dsTeK );

      // Downsample the data
      dsTrK.subsample( v );

      // Verify the synchronization
      for( unsigned int j = 0; j < dsTrK.size(); ++j )
	{
	  if( dsTrPPI.i2s( j ) != dsTrK.i2s( j ) )
	    throw std::logic_error( "The source data is not synchronized between the two views" );
	}

      // Load the external data
      if( params.alg_choice() != "joint" )
	{
	  const string fn = fnExt + sfx;
	  parseSparseFile( fn, dsExK, pfm, ',', ',', '=' );
	}

      cout<<"Space "<<i<<": Train - "<<dsTrK.size()<<" ; Test - "<<dsTeK.size()<<" ; "<<" External - "<<dsExK.size()<<" samples; ";
      cout<<nFeats( dsTrK )<<" features"<<endl;

      // Append the data to the joint set
      bool bRemoveMissing = i > 1;
      expand( dsTrV1, dsTrK, bRemoveMissing );
      expand( dsTeV1, dsTeK, bRemoveMissing );
      expand( dsExV1, dsExK, bRemoveMissing );
    }

  // Combine the views for the joint classifier
  if( params.alg_choice() == "joint" )
    {
      expand( dsTrV1, dsTrPPI, true );
      expand( dsTeV1, dsTePPI, true );
    }
  else
    {
      expand( dsTrV2, dsTrPPI, false );
      expand( dsTeV2, dsTePPI, false );
    }

  // Handle the chain classification (view 1 predictions used as view 2 features)
  if( params.alg_choice() == "chain" )
    {
      CDataSet< CSparseSample > dsTrPred, dsTePred;

      // Parse the predictions
      string sfx = "-" + boost::lexical_cast<string>( nUnlabeled ) + ".pred";
      loadV1Preds( pfxV1P, sfx, nFolds, testFold, dsTrPred, dsTePred );
      
      // Append the data to the second view
      expand( dsTrV2, dsTrPred, true );
      expand( dsTeV2, dsTePred, true );
    }

  // Load the annotations, abusing the fold parameters to load everything into dsAnnots
  const string sfxAnnots = "-annots.sdat";
  shared_ptr< CFeatMap > pfmAnnots = loadFoldSpace( fnV1, sfxAnnots, nFolds, nFolds+1, dsAnnots, dsDud );
  parseSparseFile( fnExt + sfxAnnots, dsAnnots, pfmAnnots, ',', ',', '=' );
  cout<<"Annotations data has: "<<dsAnnots.size()<<" samples, "<<nFeats( dsAnnots )<<" features"<<endl;

  // Specify the kernel and loss objects
  typedef CKernelLoss<CSparseSample> CSparseLoss;
  CDataSet<vSparseSample>::binop_t fiker = CCompositeSparseKernel(false);
  CDataSet<CSparseSample>::binop_t foker = CSparseHomKernel(true);
  CDataSet<CSparseSample>::binop_t floss = CSparseLoss(foker);
  function< double( double, double ) > fioker = CProdJointKernel();

  // Create common output-space storage
  shared_ptr< CDataSet<CSparseSample> > pO( new CDataSet<CSparseSample>( foker, floss ) );

  // Create the IO datasets
  SData data( fiker, fiker, pO, fioker );

  // Add the data samples
  if( params.alg_choice() != "view2" && params.alg_choice() != "chain" )
    {
      //      data.pTrainV1->addSets( dsTrV1, dsAnnots );
      addData( *data.pTrainV1, dsTrV1, dsAnnots, nLabeled, nUnlabeled );
      data.pTrainV1->addSets( dsExV1, dsAnnots );
      data.pTestV1->addSets( dsTeV1, dsAnnots );
    }
  if( params.alg_choice() != "view1" && params.alg_choice() != "joint" )
    {
      //      data.pTrainV2->addSets( dsTrV2, dsAnnots );
      addData( *data.pTrainV2, dsTrV2, dsAnnots, nLabeled, nUnlabeled );
      data.pTestV2->addSets( dsTeV2, dsAnnots );
    }
  data.display();

  // Datasets should not contain duplicates
  checkDuplicates( data.pTrainV1->getI() );
  checkDuplicates( data.pTrainV2->getI() );
  checkDuplicates( data.pTestV1->getI() );
  checkDuplicates( data.pTestV2->getI() );

  // Cache all kernel and loss values
  data.pTrainV1->cache();
  data.pTrainV2->cache();
  data.pTestV1->cache();
  data.pTestV2->cache();

  return data;
}

// Entry point of the code
int main( int argc, char* argv[] )
{
  // Display the usage and/or parse the parameters
  if( argc < 2 ) { cout<<"Usage: "<<argv[0]<<" <options file>"<<std::endl; return -1; }
  CStrutParams params; params.load( argv[1] );
  
  const unsigned int nFolds = 5;
  const string fnPrefix = "mview/yeast/";

  // Parse the experiment type
  bool bParamSel = false;
  if( params.exp_type() == string( "test" ) )
    cout<<"Running a test experiment"<<endl;
  else if( params.exp_type() == string( "ps" ) )
    { cout<<"Running parameter selection"<<endl; bParamSel = true; }
  else throw std::logic_error( "Unknown experiment type" );

  // Parse the fold choices
  vector< unsigned int > myFolds = params.folds();
  unsigned int testFold = myFolds[0];
  cout<<"Test fold: "<<testFold<<endl;
  int psFold = -1; 
  if( bParamSel )
    {
      psFold = myFolds[1];
      cout<<"Parameter selection test fold: "<<psFold<<endl;
    }

  if( testFold > 4 ) throw std::logic_error( "Fold index must be between 0 and 4" );
  if( psFold > 3 ) throw std::logic_error( "Parameter selection fold index must be between 0 and 3" );

  // Parse the number of sample
  vector< double > myAlgParams = params.alg_params();
  unsigned int nLabeled = static_cast<unsigned int>(myAlgParams[0]);
  unsigned int nUnlabeled = static_cast<unsigned int>(myAlgParams[1]);

  SData data = loadData( fnPrefix, nFolds, testFold, params, nLabeled, nUnlabeled );

  // Setup the prefix for saving intermediate models
  string clsfPrefix = params.log_name();
  if( *(clsfPrefix.rbegin()) != '/' ) clsfPrefix += "/";
  clsfPrefix += boost::lexical_cast<string>( testFold )
  + "-" + boost::lexical_cast<string>( nUnlabeled );
  //  string clsfPrefix = "";

  // Multiview experiment
  vector< double > vloss;
  vector< double > auroc;
  if( params.alg_choice() == "co" || params.alg_choice() == "trans" )
    {
      // Specify the parameters
      SCOSVMParams svmp;
      svmp.Cn_l = myAlgParams[2];
      svmp.Cn_u = myAlgParams[3];
      svmp.eps = 0.1;
      svmp.rmax = 10;
      svmp.fnPrefix = clsfPrefix;
      if( params.alg_choice() == "trans" )
	svmp.bTrans = true;
      else
	svmp.bTrans = false;

      // Create the COSVM
      typedef CCOSVM<vSparseSample, vSparseSample, CSparseSample> mySVM;
      shared_ptr< mySVM > clsf( new mySVM( svmp ) );

      // Training
      clsf->train( data.pTrainV1, data.pTrainV2 );
    
      // Compute the loss over the unlabeled data
      /*    double uLoss1 = 0.0;
	    double uLoss2 = 0.0;
	    //    cout<<"Infered (v1,v2) pairs: ";
	    for( unsigned int i = 0; i < nUnlabeled; i++ )
	    {
	    string s = pTrainB->i2s( i );
	    int iO = pAnnot->s2i( s );
	    std::pair<unsigned int, unsigned int> pPred = clsf->mapUnlab( s );
	    //	cout<<"("<<pPred.first<<","<<pPred.second<<")"<<endl;
	    shared_ptr< CSparseSample const > lTrue = pAnnot->getSample( iO );
	    uLoss1 += pTrainB->oloss( pPred.first, lTrue );
	    uLoss2 += pTrainI->oloss( pPred.second, lTrue );
	    }
	    //    cout<<endl;
	    cout<<"Loss / sample over unlabeled predictions in view 1: "
	    <<uLoss1/nUnlabeled<<endl;
	    cout<<"Loss / sample over unlabeled predictions in view 2: "
	    <<uLoss2/nUnlabeled<<endl;*/

      cout<<"About to the test the classifier on ";
      cout<<data.pTestV1->sizeI()<<" view1 samples and ";
      cout<<data.pTestV2->sizeI()<<" view2 samples"<<endl;

      if( data.pTestV1->getO() != data.pTestV2->getO() )
	throw std::logic_error( "Inconsistent views" );

      // Testing
      vloss = clsf->test( data.pTestV1, data.pTestV2, "", "" );

      // Compute the prediction scores
      CDataSet<CSparseSample> dsScores = predScores( *clsf, data.pTestV1->getI(), data.pTestV2->getI() );
      shared_ptr<CFeatMap const> pfmScores = dsScores.getSample(0)->getFeatMap();

      // Compute the average AUROC
      for( unsigned int fi = 0; fi < pfmScores->nFeats(); ++fi )
	{
	  string f = pfmScores->i2f( fi );

	  // Compose the (score, label) pairs
	  vector< pair< double, unsigned int > > rocData;
	  unsigned int npos = 0, nneg = 0;
	  for( unsigned int i = 0; i < dsScores.size(); ++i )
	    {
	      double score = dsScores.getSample(i)->getValue(f);

	      // Find the true label
	      string s = dsScores.i2s( i );
	      int j = data.pTestV1->getI()->s2i( s );
	      if( j < 0 ) throw std::logic_error( "Inconsistent scores dataset" );
	      if( data.pTestV2->getI()->s2i( s ) != j )
		throw std::logic_error( "Inconsistent views" );
	      unsigned int k = data.pTestV1->map(j);
	      if( data.pTestV2->map(j) != k )
		throw std::logic_error( "Inconsistent views" );
	      unsigned int lb = 0;
	      if( data.pTestV1->getO()->getSample(k)->getValue(f) != 0.0 )
		{ lb = 1; ++npos; }
	      else {++nneg;}

	      // Store the pair
	      pair< double, unsigned int > p( score, lb );
	      rocData.push_back( p );
	    }

	  if( npos == 0 || nneg == 0 ) continue;

	  // Compute the curve and the associated score
	  vector< pair< double, double > > roc = ROC( rocData );
	  auroc.push_back( AUROC( roc ) );
	}
    }
  else
    {
      // Setup the filename for saving predictions
      string fnPred = "";
      /*      if( params.alg_choice() == "view1" )
	{
	  fnPred = clsfPrefix;
	  string::size_type pos = fnPred.find( "clsf" );
	  if( pos != fnPred.npos )
	    fnPred.replace( pos, 4, "pred" );
	  fnPred += ".pred";
	  cout<<"Will save predictions to "<<fnPred<<endl;
	  }*/

      // Specify the SVM parameters
      SSSVMParams svmp;
      svmp.Cn = myAlgParams[2];
      svmp.eps = 0.1;
      svmp.nMaxQPSteps = 1000;
      svmp.fnPrefix = clsfPrefix;

      // Fetch the appropriate dataset
      shared_ptr< SData::CIOSet > pdsTrain;
      shared_ptr< SData::CIOSet > pdsTest;
      if( params.alg_choice() == "joint" || params.alg_choice() == "view1" )
	{ pdsTrain = data.pTrainV1; pdsTest = data.pTestV1; }
      else if( params.alg_choice() == "view2" || params.alg_choice() == "chain" )
	{ pdsTrain = data.pTrainV2; pdsTest = data.pTestV2; }
      else throw std::logic_error( "Unknown method type" );

      // Create and run the classifier
      CnsSSVM< vSparseSample, CSparseSample, 'm' > clsf( svmp );
      clsf.CClassifier<vSparseSample,CSparseSample>::train( pdsTrain );
      vloss = clsf.test( pdsTest, fnPred );
      
      // Compute the prediction scores
      CDataSet<CSparseSample> dsScores = predScores( clsf, pdsTest->getI() );
      shared_ptr<CFeatMap const> pfmScores = dsScores.getSample(0)->getFeatMap();
      
      // Compute the average AUROC
      for( unsigned int fi = 0; fi < pfmScores->nFeats(); ++fi )
	{
	  string f = pfmScores->i2f( fi );

	  // Compose the (score, label) pairs
	  vector< pair< double, unsigned int > > rocData;
	  unsigned int npos = 0, nneg = 0;
	  for( unsigned int i = 0; i < dsScores.size(); ++i )
	    {
	      double score = dsScores.getSample(i)->getValue(f);

	      // Find the true label
	      string s = dsScores.i2s( i );
	      int j = pdsTest->getI()->s2i( s );
	      if( j < 0 ) throw std::logic_error( "Inconsistent scores dataset" );
	      unsigned int k = pdsTest->map( j );
	      unsigned int lb = 0;
	      if( pdsTest->getO()->getSample(k)->getValue(f) != 0.0 )
		{ lb = 1; ++npos; }
	      else {++nneg;}

	      // Store the pair
	      pair< double, unsigned int > p( score, lb );
	      rocData.push_back( p );
	    }
	  
	  if( npos == 0 || nneg == 0 ) continue;

	  // Compute the curve and the associated score
	  vector< pair< double, double > > roc = ROC( rocData );
	  auroc.push_back( AUROC( roc ) );
	}
    }

  // Compute and display the average loss
  double m = std::accumulate( vloss.begin(), vloss.end(), 0.0 ) /
    static_cast<double>( vloss.size() );
  cout<<"Mean loss per sample: "<<m<<endl;

  m = std::accumulate( auroc.begin(), auroc.end(), 0.0 ) /
    static_cast<double>( auroc.size() );
  cout<<"Mean AUROC: "<<m<<endl;

  return 0;
}
