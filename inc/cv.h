// cv.h - functions that help with cross-validation partitioning
//
// by Artem Sokolov

#ifndef CV_H__INCLUDED
#define CV_H__INCLUDED

#include "types.h"

// Displays the ranges to the screen
void display( const virange_t& vTrain, const virange_t& vTest );

// Computes training and test indices given fold sizes and requested test fold
void splitCV( const vector< unsigned int >& vSizes, unsigned int iFold,
	      virange_t& vTrain, virange_t& vTest );

// Computes training and test indices given dataset size and requested test fold
// Performs equal fold split
void splitCV( const unsigned int n, const unsigned int nFolds, unsigned int iFold,
	      virange_t& vTrain, virange_t& vTest );

#endif
