// cafa-train.cpp - trains a classifier on the CAFA data
//
// by Artem Sokolov

#include "types.h"
#include "sample.h"
#include "parsers.h"
#include "io-dataset.h"
#include "clsf.h"
#include "createClsf.h"

#include <fstream>

string clsfPrefix = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/clsf/";

int main( int argc, char* argv[] )
{
  // Display usage as needed
  if( argc < 3 )
    {
      cout<<"Usage: "<<argv[0]<<" <file containing all relevant prefixes> <name of the classifier>"<<endl;
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

  pdsio->cache();

  // Create a classifier
  SSSVMParams svmp;
  svmp.Cn = 1.0;
  svmp.eps = 0.01;
  svmp.nMaxQPSteps = 1000;
  svmp.fnPrefix = clsfPrefix + sClsfName;
  shared_ptr< CClassifier<CSparseSample,CSparseSample> > psvm( new CnsSSVM<CSparseSample,CSparseSample,'m'>( svmp ) );
  psvm->train( pdsio );

  return 0;
}

