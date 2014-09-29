// cafa-train-ppi.cpp - Train a CAFA classifier based on the PPI data
//
// by Artem Sokolov

#include "types.h"
#include "sample.h"
#include "parsers.h"
#include "io-dataset.h"
#include "clsf.h"
#include "nssvm.h"

string sInputPrefix = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/final/";
string sClsfPrefix = "/s/chopin/c/proj/protfun/users/sokolov/CAFA/clsf/";

int main( int argc, char* argv[] )
{
  // Display usage
  if( argc < 3 )
    {
      cout<<"Usage: "<<argv[0]<<" <species prefix> <ontology prefix>"<<endl;
      return -1;
    }
  
  // Parse the argument and compose the filenames
  string fnPPI = sInputPrefix + argv[1] + "_ppi_train.sdat.gz";
  string fnAnnot = sInputPrefix + argv[1] + "_" + argv[2] + "_annot.sdat.gz";

  // Load the data
  shared_ptr< CDataSet< CSparseSample > > pids( new CDataSet< CSparseSample > );
  shared_ptr< CDataSet< CSparseSample > > pods( new CDataSet< CSparseSample > );
  cout<<"Parsing "<<fnPPI<<"... "; cout.flush();
  parseSparseFile( fnPPI, *pids, ',', ',', '=' );
  cout<<" parsed "<<pids->size()<<" samples, "<<nFeats( pids )<<" features"<<endl;
  cout<<"Parsing "<<fnAnnot<<"... "; cout.flush();
  parseSparseFile( fnAnnot, *pods, ',', ',', '=' );
  cout<<" parsed "<<pods->size()<<" samples, "<<nFeats( pods )<<" features"<<endl;

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
  svmp.fnPrefix = sClsfPrefix + argv[1] + "_" + argv[2];
  shared_ptr< CClassifier<CSparseSample,CSparseSample> > psvm( new CnsSSVM<CSparseSample,CSparseSample,'m'>( svmp ) );
  psvm->train( pdsio );

  return 0;
}
