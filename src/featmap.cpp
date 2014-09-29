/// \file featmap.cpp
/// \brief Implementation of CFeatMap class
/// \author Artem Sokolov

#include "featmap.h"

// Constructs a feature map out of feature ids
CFeatMap::CFeatMap( const vector< string >& fids )
  : featIDs( fids )
{
  // Build the inverse map
  for( unsigned int i = 0; i < fids.size(); i++ )
    i_featIDs[fids[i]] = i;
}

// Adds a feature to the dataset and return the new feature's index
unsigned int CFeatMap::addFeat( const string& str )
{
  // Determine if the feature already exists
  if( i_featIDs.find( str ) != i_featIDs.end() ) return i_featIDs[str];
  
  // If not, add it and update all the internal info
  unsigned int feat_i = featIDs.size();
  featIDs.push_back( str );
  i_featIDs[ str ] = feat_i;
  return feat_i;
}

// Looks up the feature name by its index
// Returns "" when index is out of bounds
const string CFeatMap::i2f( unsigned int i ) const
{
  if( i >= featIDs.size() ) return "";
  else return featIDs[i];
}

// Looks up feature's index by its name, return -1 if no such feature
int CFeatMap::f2i( const std::string& str ) const
{
  sumap_t::const_iterator iter = i_featIDs.find( str );
  if( iter == i_featIDs.end() ) return -1;
  return iter->second;
}
