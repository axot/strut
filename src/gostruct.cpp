// gostruct.cpp - prediction of GO terms using structured-outputs framework
//
// by Artem Sokolov

#include "sample.h"
#include "createClsf.h"
#include "misc.h"
#include "blast-nn.h"
#include "blast-wrap.h"
#include "params.h"
#include "parsers.h"
#include "cv.h"

#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <numeric>

#include <memory>

using std::cout;
using std::string;
using std::vector;
using std::set;
using std::shared_ptr;

string blast_hits_filename( "/s/chopin/c/proj/protfun/users/sokolov/data/foursp.blast" );

int main( int argc, char* argv[] )
{
  // Display usage and/or parse the parameters
  CGOStrutParams params;
  if( argc < 2 ) 
    { 
      cout<<"Usage: "<<argv[0]<<" <command line options or options filename>"<<std::endl; 
      params.displayHelp();
      return -1; 
    }
  params.load( argv[1] ); params.display();

  // typedefs
  typedef CKernelLoss<CSparseSample> CMyLoss;
  typedef CBLASTNN<CSparseSample, CSparseSample> CMyBNN;
  typedef CDataSet<CSparseSample> CSDataSet;

  // Load the blast hits
  string str = blast_hits_filename;
  cout<<"Loading "<<str<<"..."; cout.flush();
  shared_ptr< GO::CBLASTOutput > blast_hits( new GO::CBLASTOutput( str.c_str() ) );
  cout<<blast_hits->size()<<" entries parsed"<<std::endl;

  // Create the kernels and the loss
  CDataSet<CSparseSample>::binop_t fkeri = CSparseKernel( true );
  CDataSet<CSparseSample>::binop_t fkero = CSparseHomKernel( true );
  CDataSet<CSparseSample>::binop_t floss = CKernelLoss<CSparseSample>( fkero );
  function< double(double, double) > fioker = CProdJointKernel();

  // Create an empty input-output dataset
  typedef CIODataSet<CSparseSample, CSparseSample> CSparseIODS;
  shared_ptr< CSparseIODS > pds( new CSparseIODS( fkeri, fkero, floss, fioker ) );

  // Load the input space data
  cout<<"Loading input space... "; cout.flush();
  shared_ptr< CSDataSet > pdsi( new CSDataSet() );
  parseSparseFile( "gostruct/input.sdat", *pdsi, ',', ',', '=' );
  cout<<pdsi->size()<<" samples; "<<nFeats( pdsi )<<" features"<<std::endl;

  // Load the output space data, using a common feature map
  shared_ptr< CFeatMap > pfmap( new CFeatMap );
  vector< unsigned int > foldSizes;
  for( unsigned int i = 0; i < 4; i++ )
    {
      cout<<"Loading output space for fold "<<i<<"... "; cout.flush();
      shared_ptr< CSDataSet > pdso( new CSDataSet() );
      string filename = "gostruct/output" + boost::lexical_cast<string>(i) + ".sdat";
      parseSparseFile( filename, *pdso, pfmap, ',', ',', '=' );
      cout<<pdso->size()<<" samples; "<<nFeats( pdso )<<" features"<<std::endl;

      // DEBUG
      vector< unsigned int > v;
      for( unsigned int i = 0; i < 300; i++ ) v.push_back( i );
      pdso->subsample( v );
      // END OF DEBUG

      // Add the data to the IO dataset
      unsigned int fn = pds->addSets( pdsi, pdso );
      foldSizes.push_back( fn );
    }

  cout<<"Output space has "<<pds->getO()->size()<<" samples; "<<
    nFeats( pds->getO() )<<" features"<<endl;

  cout<<"Fold sizes: ";
  std::copy( foldSizes.begin(), foldSizes.end(),
	     std::ostream_iterator<unsigned int>( cout, " " ) );
  cout<<endl;

  // Perform the training / test split
  virange_t vTrain, vTest;
  unsigned int iFold = params.folds()[0];
  splitCV( foldSizes, iFold, vTrain, vTest );
  display( vTrain, vTest );
  shared_ptr< CSparseIODS > pTrain( static_cast<CSparseIODS*>(NULL) );
  shared_ptr< CSparseIODS > pTest( static_cast<CSparseIODS*>(NULL) );
  pds->splitTrainTest( vTrain, vTest, pTrain, pTest );
  pTrain->cache();
  pTest->cache();

  // Create the classifier
  shared_ptr< CClassifier<CSparseSample, CSparseSample> > pclsf;
  if( params.alg_choice() == string( "blast-nn" ) )
    pclsf.reset( new CMyBNN( blast_hits ) );
  else
    pclsf = createClassifier<CSparseSample, CSparseSample>( params );

  // Train the classifier
  pclsf->train( pTrain );

  // Test the classifier on the 3 OSGs
  vector< double > loss = pclsf->test( pTest );
  
  // Compute the mean loss per example
  double m = std::accumulate(loss.begin(),loss.end(),0.0) / static_cast<double>( loss.size() );

  // Report the values
  cout<<"Mean loss per test sample: "<<m<<endl;

  return 0;
}
