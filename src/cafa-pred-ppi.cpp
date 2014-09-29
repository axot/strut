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
      cout<<"<beginning index> <1 past the final index> <ontology>"<<endl;
      return -1;
    }

  // Parse the argument and compose the filenames
  string fnPPI = sInputPrefix + argv[3] + "_ppi_train.sdat.gz";
  string fnAnnot = sInputPrefix + argv[3] + "_" + argv[6] + "_annot.sdat.gz";

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

  // Load the PPI data
  shared_ptr< CFeatMap > pfmiPPI( new CFeatMap );
  shared_ptr< CDataSet< CSparseSample > > pidsPPI( new CDataSet< CSparseSample > );
  shared_ptr< CDataSet< CSparseSample > > podsPPI( new CDataSet< CSparseSample > );
  cout<<"Parsing "<<fnPPI<<"... "; cout.flush();
  parseSparseFile( fnPPI, *pidsPPI, pfmiPPI, ',', ',', '=' );
  cout<<" parsed "<<pidsPPI->size()<<" samples, "<<nFeats( pidsPPI )<<" features"<<endl;
  cout<<"Parsing "<<fnAnnot<<"... "; cout.flush();
  parseSparseFile( fnAnnot, *podsPPI, pfmo, ',', ',', '=' );
  cout<<" parsed "<<podsPPI->size()<<" samples, "<<nFeats( podsPPI )<<" features"<<endl;

  // Create the kernel/loss objects
  typedef CKernelLoss<CSparseSample> CSparseLoss;
  CDataSet<CSparseSample>::binop_t fkeri = CSparseKernel( true );
  CDataSet<CSparseSample>::binop_t fkero = CSparseHomKernel( true );
  CDataSet<CSparseSample>::binop_t floss = CSparseLoss( fkero );
  function< double( double, double ) > fkerio=  CProdJointKernel();

  // Create an IO mapping
  typedef CIODataSet< CSparseSample, CSparseSample > CIOSet;
  shared_ptr< CIOSet > pdsio( new CIOSet( fkeri, fkero, floss, fkerio ) );
  pdsio->addSets( pids, pods );
  shared_ptr< CIOSet > pdsioPPI( new CIOSet( fkeri, fkero, floss, fkerio ) );
  pdsioPPI->addSets( pidsPPI, podsPPI );

  cout<<"Cross-species";
  cout<<"  Input space has "<<pdsio->sizeI()<<" samples, ";
  cout<<nFeats( pdsio->getI() )<<" features"<<endl;
  cout<<"  Output space has "<<pdsio->sizeO()<<" samples, ";
  cout<<nFeats( pdsio->getO() )<<" features"<<endl;

  cout<<"PPI";
  cout<<"  Input space has "<<pdsioPPI->sizeI()<<" samples, ";
  cout<<nFeats( pdsioPPI->getI() )<<" features"<<endl;
  cout<<"  Output space has "<<pdsioPPI->sizeO()<<" samples, ";
  cout<<nFeats( pdsioPPI->getO() )<<" features"<<endl;

  pdsio->cache();
  pdsioPPI->cache();

  // Create a classifier
  SSSVMParams svmp;
  svmp.Cn = 1.0;
  svmp.eps = 0.01;
  svmp.nMaxQPSteps = 1000;
  svmp.fnPrefix = clsfPrefix + sClsfName;
  shared_ptr< CnsSSVM<CSparseSample,CSparseSample,'m'> > psvm( new CnsSSVM<CSparseSample,CSparseSample,'m'>( svmp ) );

  SSSVMParams svmpPPI;
  svmpPPI.Cn = 1.0;
  svmpPPI.eps = 0.01;
  svmpPPI.nMaxQPSteps = 1000;
  svmpPPI.fnPrefix = clsfPrefix + argv[3] + "_" + argv[6];
  shared_ptr< CnsSSVM<CSparseSample,CSparseSample,'m'> > psvmPPI( new CnsSSVM<CSparseSample,CSparseSample,'m'>( svmpPPI ) );

  // Load the latest model
  int curIter = psvm->preload( pdsio );
  cout<<"Preloaded iteration "<<curIter<<endl;

  curIter = psvmPPI->preload( pdsioPPI );
  cout<<"PPI: Preloaded iteration "<<curIter<<endl;

  // Load the test data
  string fnTest = sInputPrefix + argv[3] + "_test.sdat.gz";
  cout<<"Loading test data... "; cout.flush();
  shared_ptr< CDataSet< CSparseSample > > pTest( new CDataSet< CSparseSample >(fkeri) );
  parseSparseFile( fnTest, *pTest, pfmi, ',', ',', '=' );
  cout<<"parsed "<<pTest->size()<<" samples, "<<nFeats( pTest )<<" features"<<endl;

  // Load the PPI test data
  string fnTestPPI = sInputPrefix + argv[3] + "_ppi_test.sdat.gz";
  cout<<"Loading PPI test data... "; cout.flush();
  shared_ptr< CDataSet< CSparseSample > > pTestPPI( new CDataSet< CSparseSample >(fkeri) );
  parseSparseFile( fnTestPPI, *pTestPPI, pfmiPPI, ',', ',', '=' );
  cout<<"parsed "<<pTestPPI->size()<<" samples, "<<nFeats( pTestPPI )<<" features"<<endl;

  // Limit the scope to the examples of interest
  unsigned int iBegin = boost::lexical_cast< unsigned int >( argv[4] );
  unsigned int iEnd = boost::lexical_cast< unsigned int >( argv[5] );
  vector< unsigned int > v;
  for( unsigned int i = iBegin; i < iEnd; ++i )
    v.push_back( i );
  pTestPPI->subsample( v );
  cout<<"Classifying ["<<iBegin<<","<<iEnd<<")"<<endl;

  // Perform the predictions
  string fnPred = sPredPrefix + argv[3] + "_" + argv[6] + "_" + argv[4] + "-" + argv[5];
  cout<<"Saving predictions to "<<fnPred<<endl;
  std::ofstream ofs( fnPred.c_str() );

  // Construct a map from PPI output space to cross-species output space
  uumap_t annotMap;
  for( unsigned int i = 0; i < pdsioPPI->sizeO(); ++i )
    {
      // Find the sample
      shared_ptr< CSparseSample const > p = pdsioPPI->getO()->getSample( i );
      int j = pdsio->getO()->findSample( p );
      if( j < 0 ) 
	throw std::logic_error( "Unable to sync output spaces" );
      
      // Store the entry
      annotMap[i] = j;
    }

  pdsio->cacheIExternal( pTest );
  pdsioPPI->cacheIExternal( pTestPPI );
  shared_ptr< const CDataSet<CSparseSample> > pO1 = pdsioPPI->getO();
  shared_ptr< const CDataSet<CSparseSample> > pO2 = pdsio->getO();
  for( unsigned int i = 0; i < pTestPPI->size(); ++i )
    {
      // Display progress
      if( i % 100 == 0 ) {cout<<"."; cout.flush();}

      string s = pTestPPI->i2s( i );
      int ii = pTest->s2i( s );
      if( ii < 0 ) throw std::logic_error( "Sample doesn't exist in the cross-species space" );

      // Compose the sample pairs
      std::pair< const CDataSet<CSparseSample>&, unsigned int > smpl1( *pTestPPI, i );
      std::pair< const CDataSet<CSparseSample>&, unsigned int > smpl2( *pTest, ii );

      // Traverse the output examples
      unsigned int amax = 0;
      double max = -1.0 * std::numeric_limits<double>::infinity();
      for( unsigned int j = 0; j < pO1->size(); ++j )
	{
	  double v1 = psvmPPI->f( smpl1, j );
	  double v2 = psvm->f( smpl2, annotMap[j] );
	  double v = 0.2 * v1 + 0.8 * v2;
	  if( v > max ) { max = v; amax = j; }
	}

      // Save the prediction
      ofs<<pTestPPI->i2s(i)<<","<<*(pO1->getSample( amax ));
    }
  cout<<endl;
  ofs.close();

  return 0;
}
