// cv.cpp - functions that help with cross-validation partitioning
//
// by Artem Sokolov

#include "cv.h"

// Displays the ranges to the screen
void display( const virange_t& vTrain, const virange_t& vTest )
{
  cout<<"Training ranges :";
  for( unsigned int i = 0; i < vTrain.size(); i++ )
    cout<<" ["<<vTrain[i].first<<" "<<vTrain[i].second<<")";
  cout<<endl;
  cout<<"Test ranges :";
  for( unsigned int i = 0; i < vTest.size(); i++ )
    cout<<" ["<<vTest[i].first<<" "<<vTest[i].second<<")";
  cout<<endl;
}

// Computes training and test indices given fold sizes and requested test fold
void splitCV( const vector< unsigned int >& vSizes, unsigned int iFold,
	      virange_t& vTrain, virange_t& vTest )
{
  // Verify the arguments
  if( iFold >= vSizes.size() ) throw std::logic_error( "Fold index is out of range" );

  // Prep the output
  vTrain.clear();
  vTest.clear();

  // Traverse the folds
  unsigned int iLeft = 0;
  unsigned int iRight = 0;
  for( unsigned int i = 0; i < vSizes.size(); i++ )
    {
      // Compose the range
      iRight = iLeft + vSizes[i];
      irange_t r( iLeft, iRight );

      // Handle the test fold
      if( i == iFold )
	vTest.push_back( r );
      else
	vTrain.push_back( r );

      // Update the bound
      iLeft = iRight;
    }
}

// Computes training and test indices given dataset size and requested test fold
void splitCV( const unsigned int n, const unsigned int nFolds, unsigned int iFold,
	      virange_t& vTrain, virange_t& vTest )
{
  // Split the n examples into equal folds
  unsigned int nWhole = n / nFolds;
  unsigned int nRemain = n - nWhole * nFolds;
  
  // Compose the fold sizes
  vector< unsigned int > vSizes;
  for( unsigned int i = 0; i < nFolds; i++ )
    {
      unsigned int ni = nWhole;
      if( i < nRemain ) ni++;
      vSizes.push_back( ni );
    }

  splitCV( vSizes, iFold, vTrain, vTest );
}
