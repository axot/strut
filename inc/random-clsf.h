// -*-c++-*-
// \file random-clsf.h
// \brief Random classifier for structured output prediction
// \author Artem Sokolov

#ifndef RANDOM_CLSF_H__INCLUDED
#define RANDOM_CLSF_H__INCLUDED

#include "misc.h"

#include <memory>

using std::shared_ptr;

/// Random classifier for structured output prediction
template< typename _I, typename _O >
class CRandomClassifier : public CClassifier<_I, _O>
{
public:
  /// Constructor
  explicit CRandomClassifier()
    : CClassifier<_I, _O>( "Random" ) {}

private:
  /// Trains the classifier  
  virtual void train() {}

public:
  /// Clear the learned information for the classifier
  virtual void clear() {}

  /// Generic implementation of the compatibility function is not used
  virtual double f( pair< const CDataSet<_I>&, unsigned int > x, unsigned int dso ) const;
};

template< typename _I, typename _O >
double CRandomClassifier<_I,_O>::f( pair< const CDataSet<_I>&, unsigned int > x, unsigned int yj ) const
{
  // Generate a random value from (0, 1)
  boost::uniform_real<> unif_range( 0, 1 );
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rng(g_RNG, unif_range);
  return rng();
}

#endif

