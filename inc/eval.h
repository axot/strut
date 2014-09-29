// -*-c++-*-
/// \file eval.h
/// \brief Miscellaneous functions for evaluation
/// \author Artem Sokolov

#ifndef EVAL_H__INCLUDED
#define EVAL_H__INCLUDED

#include "types.h"

#include <map>

using std::multimap;

/// Given a ranking map of the form (confidence, truth) returns threshold associated with the best balanced rate
double bestSRate( const multimap<double, bool>& m );

/// Computes the ROC curve for a collection of (score, label) pairs
/** Returns a collection of (FP, TP) points that define the ROC curve
 */
vector< pair< double, double > > ROC( const vector< pair< double, unsigned int > >& data );

/// Computes the area under ROC for a vector of scores
double AUROC( const vector< pair<double, double> >& curve );


#endif
