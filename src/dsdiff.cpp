// dsdiff.cpp - Determines if two datasets are different
//
// by Artem Sokolov

#include "sample.h"
#include "parsers.h"

#include <iostream>
#include <memory>

using std::cout;
using std::endl;
using std::shared_ptr;

int main( int argc, char* argv[] )
{
  // Argument verification
  if( argc < 3 )
    {
      cout<<"Usage: "<<argv[0]<<" dataset1 dataset2"<<endl;
      return -1;
    }

  // Load the two datasets
  shared_ptr< CDataSet< CSparseSample > > pds1( new CDataSet< CSparseSample > );
  parseSparseFile( argv[1], *pds1, ',', ' ', ':' );
  cout<<"Dataset "<<argv[1]<<" has "<<pds1->size()<<" samples; "<<nFeats(pds1)<<" features"<<endl;
  shared_ptr< CDataSet< CSparseSample > > pds2( new CDataSet< CSparseSample > );
  parseSparseFile( argv[2], *pds2, ',', ' ', ':' );
  cout<<"Dataset "<<argv[2]<<" has "<<pds2->size()<<" samples; "<<nFeats(pds2)<<" features"<<endl;

  bool result = true;

  // Compare the number of samples and features
  if( pds1->size() != pds2->size() ) result = false;
  if( nFeats( pds1 ) != nFeats( pds2 ) ) result = false;

  // Traverse the samples
  unsigned int n = pds1->size();
  for( unsigned int i = 0; i < n; i++ )
    {
      // Display progress
      if( (i % 100) == 0 ) {cout<<"."; cout.flush();}
      
      // Find the sample in the other dataset
      string s = pds1->i2s( i );
      int i2 = pds2->s2i( s );
      if( i2 < 0 ) {result = false; break; }

      // Fetch the two samples and their feature maps
      shared_ptr< const CSparseSample > p1 = pds1->getSample( i );
      shared_ptr< const CFeatMap > pfm1 = p1->getFeatMap();
      shared_ptr< const CSparseSample > p2 = pds2->getSample( i2 );
      shared_ptr< const CFeatMap > pfm2 = p2->getFeatMap();

      if( pfm1->nFeats() != pfm2->nFeats() ) { result = false; break; }

      // Traverse the features
      unsigned int m = pfm1->nFeats();
      for( unsigned int j = 0; j < m; j++ )
	{
	  // Find the feature in the other map
	  string f = pfm1->i2f( j );
	  int j2 = pfm2->f2i( f );
	  if( j2 < 0 ) {result = false; break; }

	  // Compare the values
	  if( p1->getValue( j ) != p2->getValue( j2 ) )
	    { result = false; break; }
	}

      if( result == false ) break;
    }
  cout<<endl;

  if( result == true )
    {
      cout<<"Datasets "<<argv[1]<<" and "<<argv[2]<<" are identical"<<endl;
      return 0;
    }
  else
    {
      cout<<"Datasets "<<argv[1]<<" and "<<argv[2]<<" differ"<<endl;
    }
}
