// mousefunc.cpp - experiments with the mousefunc dataset
//
// by Artem Sokolov

#include "sample.h"
#include "params.h"
#include "createClsf.h"
#include "parsers.h"
#include "io-dataset.h"
#include "cv.h"

#include <iostream>
#include <fstream>
#include <set>
#include <numeric>

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
#include <boost/spirit/include/classic_increment_actor.hpp>
#include <boost/spirit/include/classic_clear_actor.hpp>

#include <boost/program_options.hpp>

namespace spirit = boost::spirit::classic;
namespace po = boost::program_options;

using std::cout;
using std::endl;
using std::set;

// Loads a new feature space and appends it to the multiset
void loadFeatureSpace( string filename, CDataSet< vSparseSample >& ds )
{
  // Load the data
  cout<<"Loading "<<filename<<"... "; cout.flush();
  CDataSet<CSparseSample> d;
  parseSparseFile( filename, d, ',', ',', '=' );
  cout<<"Loaded "<<d.size()<<" samples; "<<nFeats( d )<<" features"<<endl;

  // Append the Zhang dataset to the multiset
  expand( ds, d );
  cout<<"  The combined dataset is now at "<<ds.size()<<" samples"<<endl;
}

// Loads a particular ontology for the output space
void loadOntologySpace( const string& sFull, const string& sAbbrev,
			CDataSet< CSparseSample >& dsoTr, CDataSet< CSparseSample >& dsoTe,
			shared_ptr< CFeatMap > pfm )
{
  CDataSet< CSparseSample > d1;
  CDataSet< CSparseSample > d2;

  string fnTrain = "mousefunc/" + sAbbrev + "_train_c.tdel";
  string fnTest = "mousefunc/" + sAbbrev + "_test_c.tdel";
  
  // Load the training data
  cout<<"Loading "<<sFull<<"... "<<endl;
  parseTabDelFile( fnTrain, d1 );
  expand( dsoTr, pfm, d1 );

  // Load the test data
  cout<<"Loading "<<sFull<<"(test)... "<<endl;
  parseTabDelFile( fnTest, d2 );
  expand( dsoTe, pfm, d2 );
}

int main( int argc, char* argv[] )
{
  // Display usage and/or parse the parameters
  if( argc < 2 ) { cout<<"Usage: "<<argv[0]<<" <options file>"<<endl; return -1; }
  CGOStrutParams params; params.load( argv[1] ); params.display();

  typedef CKernelLoss<CSparseSample> CSparseKLoss;
  CDataSet<vSparseSample>::binop_t fcker = CCompositeSparseKernel( false );
  CDataSet<CSparseSample>::binop_t fker = CSparseKernel( true );
  CDataSet<CSparseSample>::binop_t floss = CSparseKLoss( fker );
  function< double( double, double ) > fioker = CProdJointKernel();

  ///////////////////////////////

  CDataSet< vSparseSample > dsi;

  // Load the datasets
  loadFeatureSpace( "mousefunc/ge_zhang.sdat.gz", dsi );
  loadFeatureSpace( "mousefunc/ge_su.sdat.gz", dsi );
  loadFeatureSpace( "mousefunc/i_adj.sdat", dsi );
  loadFeatureSpace( "mousefunc/dd_pfam.sdat", dsi );
  loadFeatureSpace( "mousefunc/dd_inter.sdat", dsi );
  loadFeatureSpace( "mousefunc/phylo.sdat", dsi );

  shared_ptr< CFeatMap > pfm( new CFeatMap );
  CDataSet<CSparseSample> dsoTr;
  CDataSet<CSparseSample> dsoTe;

  // Load the appropriate set of labels
  if( GO::hasMF( params.ontology() ) )
    loadOntologySpace( "molecular function", "mf", dsoTr, dsoTe, pfm );
  if( GO::hasBP( params.ontology() ) )
    loadOntologySpace( "biological process", "bp", dsoTr, dsoTe, pfm );
  if( GO::hasCC( params.ontology() ) )
    loadOntologySpace( "cellular compoennt", "cc", dsoTr, dsoTe, pfm );

  if( dsoTr.size() < 1 || dsoTe.size() < 1 )
    throw std::runtime_error( "No labels loaded" );

  // Append hierarchy to the samples
  cout<<"Training data output space has "<<dsoTr.size()<<" samples; "<<nFeats( dsoTr )<<" features"<<endl;
  cout<<"Test     data output space has "<<dsoTe.size()<<" samples; "<<nFeats( dsoTe )<<" features"<<endl;

  // DEBUG
  /*  vector< unsigned int > v;
  for( unsigned int i = 0; i < 400; i++ )
    v.push_back( i );
    dsoTr.subsample( v );*/
  // END OF DEBUG

  // Combine the data to form a joint output space
  typedef CIODataSet<vSparseSample, CSparseSample> CIOSet;
  shared_ptr< CIOSet > pdsio( new CIOSet(fcker, fker, floss, fioker) );
  pdsio->addSets( dsi, dsoTr );
  pdsio->addSets( dsi, dsoTe );

  cout<<"Input space has "<<pdsio->sizeI()<<" samples"<<endl;
  cout<<"Output space has "<<pdsio->sizeO()<<" samples, "<<nFeats( pdsio->getO() )<<" features"<<endl;

  string ont;
  if( GO::hasMF( params.ontology() ) ) ont += "mf";
  if( GO::hasBP( params.ontology() ) ) ont += "bp";
  if( GO::hasCC( params.ontology() ) ) ont += "cc";
  string fnScores = "mfunc-score-" + params.alg_choice() + "-" + ont + ".sdat";
  cout<<"Saving prediction scores to "<<fnScores<<endl;

  // Re-split the data into training and test sets, both using the same output space now

  // Retrieve the fold sizes
  vector< unsigned int > foldSizes;
  foldSizes.push_back( pdsio->sizeI() - dsoTe.size() );
  foldSizes.push_back( dsoTe.size() );
  
  // Create the ranges
  virange_t vTrain, vTest;
  splitCV( foldSizes, 1, vTrain, vTest );
  display( vTrain, vTest );
  
  // Perform the split
  shared_ptr< CIOSet > pTrain;
  shared_ptr< CIOSet > pTest;
  cout<<"About to split"<<endl;
  pdsio->splitTrainTest( vTrain, vTest, pTrain, pTest );
  cout<<"Finished the split"<<endl;
  pTrain->cache();
  pTest->cache();

  // Create and add the appropriate test classifier
  shared_ptr< CClassifier<vSparseSample, CSparseSample> > pclsf = createClassifier<vSparseSample, CSparseSample>( params );

  // Train and test the classifier
  pclsf->train( pTrain );
  vector< double > loss = pclsf->test( pTest );

  // Report the results
  double m = std::accumulate( loss.begin(), loss.end(), 0.0 ) /
    static_cast<double>( loss.size() );
  cout<<"Mean loss per test sample: "<<m<<endl;

  pTrain->cacheIExternal( pTest->getI() );
  predScores( *pclsf, pTest->getI(), fnScores );

  return 0;
}
