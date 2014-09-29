// -*-c++-*-
/// \file featmap.h
/// \brief Interface to the CFeatMap class
/// \author Artem Sokolov

#ifndef FEATMAP_H__INCLUDED
#define FEATMAP_H__INCLUDED

#include "types.h"

class CFeatMap
{
public:
  /// Constructor
  CFeatMap() {}

  /// Constructor
  CFeatMap( const vector< string >& fids );

  /// Destructor
  ~CFeatMap() {}

private:
  /// Maps an index to feature ID
  vector< string > featIDs;

  /// Maps a feature ID to its index
  sumap_t i_featIDs;

public:
  /// Adds a feature to the dataset and return the new feature's index
  /** If the feature ID already exists, the function returns its index
   */
  unsigned int addFeat( const string& str );

  /// Looks up the name of a feature by its index
  /** Returns "" when index is out of bounds
   */
  const string i2f( unsigned int i ) const;

  /// Looks up feature's index by its name, return -1 if no such feature
  int f2i( const string& str ) const;

  /// Returns all feature IDs
  const vector< string > getFeatureIDs() const { return featIDs; }

  /// Returns the total number of features
  unsigned int nFeats() const { return featIDs.size(); }
};

#endif
