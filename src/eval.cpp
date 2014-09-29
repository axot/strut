/// \file eval.cpp
/// \brief Implementation of various evaluation functions
/// \author Artem Sokolov

#include "eval.h"

using std::pair;

// Given a ranking map of the form (confidence, truth) returns threshold associated with the best balanced rate
double bestSRate( const multimap<double, bool>& m )
{
  // Count the total number of positives and negatives
  unsigned int nPos = 0, nNeg = 0;
  for( multimap<double, bool>::const_iterator iter = m.begin(); iter != m.end(); iter++ )
    {
      if( iter->second == true ) nPos++; else nNeg++;
    }

  // Start off with all correct positives and no correct negatives
  unsigned int nCorPos = nPos, nCorNeg = 0;

  // Traverse the map and figure out the best thershold
  double best = 0.49; double thresh = 0.0;
  for( multimap<double, bool>::const_iterator iter = m.begin(); iter != m.end(); iter++ )
    {
      if( iter->second == true ) nCorPos--;
      else nCorNeg++;

      // Account for equality with the next value
      multimap<double,bool>::const_iterator iter_new = iter; iter_new++;
      if( iter_new != m.end() && iter_new->first == iter->first ) continue;

      // Compute the balanced success rate and compare it to the best seen so far
      double fPos = static_cast<double>(nCorPos) / static_cast<double>(nPos);
      double fNeg = static_cast<double>(nCorNeg) / static_cast<double>(nNeg);
      double val = 0.5 * fPos + 0.5 * fNeg;

      // Compare to the best-so-far
      if( val > best ) { best = val; thresh = iter->first; }
    }

  return thresh;
}

/// Computes the ROC curve for a collection of (score, label) pairs
/** Returns a collection of (FP, TP) points that define the ROC curve
 */
vector< pair< double, double > > ROC( const vector< pair< double, unsigned int > >& data )
{
  multimap< double, unsigned int > mymap;

  // Stick (score, label) pairs into a map and compute the number of positive and negative labels
  unsigned int npos = 0, nneg = 0;
  for( unsigned int i = 0; i < data.size(); i++ )
    {
      if( data[i].second == 0 ) ++nneg;
      else ++npos;
      mymap.insert( data[i] );
    }

  // Traverse the map and compute the ROC curve
  vector< pair<double, double> > curve;
  multimap<double, unsigned int>::reverse_iterator iter = mymap.rbegin();
  unsigned int rpos = 0, rneg = 0;	// Running counts
  curve.push_back( pair<double,double>(0.0, 0.0) );
  while( iter != mymap.rend() )
    {
      // Find all the labels associated with this score
      double curVal = iter->first;
      vector< unsigned int > curLabels; curLabels.push_back( iter->second );
      ++iter;
      while( iter != mymap.rend() && iter->first == curVal )
	{
	  curLabels.push_back( iter->second );
	  ++iter;
	}

      // Update the running counts
      for( unsigned int i = 0; i < curLabels.size(); i++ )
	{
	  if( curLabels[i] == 0 ) rneg++;
	  else rpos++;
	}

      // Add the new point to the ROC curve
      double fp = static_cast<double>( rneg ) / static_cast<double>( nneg );
      double tp = static_cast<double>( rpos ) / static_cast<double>( npos );
      pair<double, double> pt( fp, tp );
      curve.push_back( pt );
    }
  
  return curve;
}

// Computes the area under ROC for a vector of scores
double AUROC( const vector< pair< double, double > >& curve )
{
  // Compute the area under the curve using the trapezoidal rule
  double res = 0.0;
  for( unsigned int i = 1; i < curve.size(); i++ )
    {
      double x1 = curve[i-1].first; double y1 = curve[i-1].second;
      double x2 = curve[i].first;   double y2 = curve[i].second;
      double dx = x2 - x1;
      if( dx <= 0 ) continue;
      res += 0.5 * dx * (y1 + y2);
    }

  return res;
}



