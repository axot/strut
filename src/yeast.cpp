// yeast.cpp - experiments with the yeast organisms
//
// by Artem Sokolov

#include "sample.h"
#include "cv.h"
#include "parsers.h"
#include "createClsf.h"

#include <numeric>

const string fnPrefix = "yeast/mf_";

int main( int argc, char* argv[] )
{
  // Display usage and/or parse the parameters
  if( argc < 2 ) { cout<<"Usage: "<<argv[0]<<" <options file>"<<endl; return -1; }
  CStrutParams params; params.load( argv[1] ); params.display();

  // Create the kernel/loss objects
  typedef CDataSet<CSparseSample> CSparseDataSet;
  typedef CKernelLoss<CSparseSample> CSparseLoss;
  CDataSet<vSparseSample>::binop_t fkeri = CCompositeSparseKernel( false );
  CDataSet<CSparseSample>::binop_t fkero = CSparseHomKernel( true );
  CDataSet<CSparseSample>::binop_t floss = CSparseLoss( fkero );
  function< double( double, double ) > fkerio = CProdJointKernel();

  // Load the input-space data
  vector< string > vfn;
  vfn.push_back( fnPrefix + "blast.sdat" );
  vfn.push_back( fnPrefix + "ppi.sdat" );
  vfn.push_back( fnPrefix + "local.sdat" );
  vfn.push_back( fnPrefix + "transmem.sdat" );
  CDataSet< vSparseSample > dsi = loadKernels( vfn );  

  // Load the output-space data
  CSparseDataSet dso;
  string fnAnnot = fnPrefix + "annot.sdat";
  parseSparseFile( fnAnnot, dso, ',', ',', '=' );
  cout<<"Loaded "<<dso.size()<<" annotations"<<endl;

  // DEBUG
  vector< unsigned int > v;
  for( unsigned int i = 0; i < 300; ++i ) v.push_back( i );
  dso.subsample( v );
  // END OF DEBUG

  // Create an IO-dataset and randomize the samples
  typedef CIODataSet<vSparseSample, CSparseSample> CIOSet;
  shared_ptr< CIOSet > pdsio( new CIOSet( fkeri, fkero, floss, fkerio ) );
  pdsio->addSets( dsi, dso );
  pdsio->random_shuffle();

  // Perform the fold split
  const unsigned int nFolds = 5;
  unsigned int testFold = params.folds()[0];
  virange_t vTrain, vTest;
  splitCV( dso.size(), nFolds, testFold, vTrain, vTest );
  display( vTrain, vTest );
  
  // Split the dataset
  shared_ptr< CIOSet > pTrain, pTest;
  pdsio->splitTrainTest( vTrain, vTest, pTrain, pTest );

  cout<<"Training data has "<<pTrain->sizeI()<<" samples"<<endl;
  cout<<"Test data has "<<pTest->sizeI()<<" samples"<<endl;
  cout<<"Output space has "<<pTrain->sizeO()<<" samples; "<<nFeats( pTrain->getO() )<<" features"<<endl;

  cout<<"---"<<endl;
  shared_ptr< CSparseDataSet > pO = std::const_pointer_cast< CSparseDataSet >( pTrain->getO() );
  shared_ptr< CSparseSample > p = pO->getSampleMod( 1 );
  shared_ptr< CSparseSample > q = pO->getSampleMod( 2 );

  cout<<"p - "<<*p;
  cout<<"q - "<<*q;

  (*p) += (*q);

  cout<<"p + q - "<<*p;
  
  return 0;

  // Cache the kernel values
  pTrain->cache();
  pTest->cache();

  // Create the classifier and apply it to the data
  SSSVMParams svmp;
  svmp.Cn = params.alg_params()[0];
  svmp.eps = 0.01;
  svmp.nMaxQPSteps = 1000;
  svmp.fnPrefix = params.log_name();
  shared_ptr< CClassifier<vSparseSample,CSparseSample> > pclsf( new CnsSSVM<vSparseSample,CSparseSample,'m'>( svmp ) );
  pclsf->train( pTrain );
  vector< double > loss = pclsf->test( pTest, "yeast-debug.pred" );

  // Compute the mean loss per example
  double m = std::accumulate( loss.begin(), loss.end(), 0.0 ) / static_cast<double>( loss.size() );
  cout<<"Mean loss per test sample: "<<m<<endl;

  return 0;
}
