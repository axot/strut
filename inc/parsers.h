// -*-c++-*-
/// \file parsers.h
/// \brief headers for a collection of parsers
/// \author Artem Sokolov

#ifndef PARSERS_H__INCLUDED
#define PARSERS_H__INCLUDED

#include "types.h"
#include "dataset.h"
#include "featmap.h"
#include "sample.h"

// Parses a file in sparse format
shared_ptr< CFeatMap > parseSparseFile( const string& filename,
					CDataSet<CSparseSample>& ds,
					char cSIDSep, char cPairSep, char cFVSep );

// Parses a file in sparse format and uses the provided feature map
void parseSparseFile( const string& filename,
		      CDataSet<CSparseSample>& ds,
		      shared_ptr< CFeatMap > pfmap,
		      char cSIDSep,  // separates sample ID from the rest
		      char cPairSep, // separates feat-val pairs
		      char cFVSep    // separates feature ID from value
		      );

// Parses a file in tab-delimited format
void parseTabDelFile( const string& filename,
		      CDataSet< CSparseSample >& ds );

/// Loads a collection of files, one per kernel space
CDataSet<vSparseSample> loadKernels( vector< string >& filenames );

#endif
