// cafa-pred.cpp - Generates predictions on the CAFA data
//
// by Artem Sokolov

#include "types.h"
#include "sample.h"
#include "parsers.h"
#include "io-dataset.h"
#include "clsf.h"
#include "createClsf.h"

#include <fstream>

string sInputPrefix = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/final/";
string clsfPrefix = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/clsf/";
string sPredPrefix = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/pred/";

int main( int argc, char* argv[] )
{
  // Display usage as needed
  if( argc < 6 )
    {
      cout<<"Usage: "<<argv[0]<<" <file containing all relevant prefixes> <name of the classifier> <test species prefix> ";
      cout<<"<beginning index> <1 past the final index>"<<endl;
      return -1;
    }

  string sClsfName( argv[2] );
  cout<<"Classifier name: "<<sClsfName<<endl;

  // Load the prefixes
  string s;
  vector< string > vPrefixes;
  std::ifstream ifs( argv[1] );
  for( std::getline( ifs, s ); ifs.fail() == false; std::getline( ifs, s ) )
    vPrefixes.push_back( s );

  cout<<"Using the following prefixes: "<<endl;
  for( unsigned int i = 0; i < vPrefixes.size(); ++i )
    cout<<vPrefixes[i]<<endl;

  // Load the data
  shared_ptr< CFeatMap > pfmi( new CFeatMap );
  shared_ptr< CFeatMap > pfmo( new CFeatMap );
  shared_ptr< CDataSet< CSparseSample > > pids( new CDataSet< CSparseSample > );
  shared_ptr< CDataSet< CSparseSample > > pods( new CDataSet< CSparseSample > );
  for( unsigned int i = 0; i < vPrefixes.size(); ++i )
    {
      // Compose the filename
      string ifn = vPrefixes[i] + "_train.sdat.gz";
      string ofn = vPrefixes[i] + "_annot.sdat.gz";

      // Parse the input file
      cout<<"Parsing "<<ifn<<"... "; cout.flush();
      parseSparseFile( ifn, *pids, pfmi, ',', ',', '=' );
      cout<<" input-space now has "<<pids->size()<<" samples, ";
      cout<<nFeats( pids )<<" features"<<endl;

      // Parse the output file
      cout<<"Parsing "<<ofn<<"... "; cout.flush();
      parseSparseFile( ofn, *pods, pfmo, ',', ',', '=' );
      cout<<" output-space now has "<<pods->size()<<" samples, ";
      cout<<nFeats( pods )<<" features"<<endl;
    }

  // Create the kernel/loss objects
  typedef CKernelLoss<CSparseSample> CSparseLoss;
  CDataSet<CSparseSample>::binop_t fkeri = CSparseKernel( true );
  CDataSet<CSparseSample>::binop_t fkero = CSparseHomKernel( true );
  CDataSet<CSparseSample>::binop_t floss = CSparseLoss( fkero );
  function< double( double, double ) > fkerio = CProdJointKernel();

  // Create an IO mapping
  typedef CIODataSet< CSparseSample, CSparseSample > CIOSet;
  shared_ptr< CIOSet > pdsio( new CIOSet( fkeri, fkero, floss, fkerio ) );
  pdsio->addSets( pids, pods );

  cout<<"Input space has "<<pdsio->sizeI()<<" samples, ";
  cout<<nFeats( pdsio->getI() )<<" features"<<endl;
  cout<<"Output space has "<<pdsio->sizeO()<<" samples, ";
  cout<<nFeats( pdsio->getO() )<<" features"<<endl;

  // Create a classifier
  SSSVMParams svmp;
  svmp.Cn = 1.0;
  svmp.eps = 0.01;
  svmp.nMaxQPSteps = 1000;
  svmp.fnPrefix = clsfPrefix + sClsfName;
  shared_ptr< CnsSSVM<CSparseSample,CSparseSample,'m'> > psvm( new CnsSSVM<CSparseSample,CSparseSample,'m'>( svmp ) );

  // Load the latest model
  int curIter = psvm->preload( pdsio );
  cout<<"Preloaded iteration "<<curIter<<endl;

  // Load the test data
  string fnTest = sInputPrefix + argv[3] + "_test.sdat.gz";
  cout<<"Loading test data... "; cout.flush();
  shared_ptr< CDataSet< CSparseSample > > pTest( new CDataSet< CSparseSample >(fkeri) );
  parseSparseFile( fnTest, *pTest, pfmi, ',', ',', '=' );
  cout<<"parsed "<<pTest->size()<<" samples, "<<nFeats( pTest )<<" features"<<endl;

  // Limit the scope to the examples of interest
  unsigned int iBegin = boost::lexical_cast< unsigned int >( argv[4] );
  unsigned int iEnd = boost::lexical_cast< unsigned int >( argv[5] );
  vector< unsigned int > v;
  for( unsigned int i = iBegin; i < iEnd; ++i )
    v.push_back( i );
  pTest->subsample( v );
  cout<<"Classifying ["<<iBegin<<","<<iEnd<<")"<<endl;

  // Perform the predictions
  string fnPred = sPredPrefix + argv[2] + "_" + argv[3] + "_" + argv[4] + "-" + argv[5];
  //  cout<<"Saving predictions to "<<fnPred<<endl;
  //  std::ofstream ofs( fnPred.c_str() );
  cout<<"Classifying... "<<endl;
  psvm->predict( pTest, fnPred );

  /*  
  pdsio->cacheIExternal( pTest );
  shared_ptr< const CDataSet<CSparseSample> > pO = pdsio->getO();
  for( unsigned int i = 0; i < pTest->size(); ++i )
    {
      // Display progress
      if( i % 100 == 0 ) {cout<<"."; cout.flush();}

      // Compose the sample pair
      std::pair< const CDataSet<CSparseSample>&, unsigned int > smpl( *pTest, i );

      // Traverse the output examples
      unsigned int amax = 0;
      double max = psvm->f( smpl, 0 );
      for( unsigned int j = 1; j < pO->size(); ++j )
	{
	  double v = psvm->f( smpl, j );
	  if( v > max ) { max = v; amax = j; }
	}

      // Save the prediction
      ofs<<pTest->i2s(i)<<","<<*(pO->getSample( amax ));
    }
  cout<<endl;
  ofs.close();
  */

  return 0;
}
