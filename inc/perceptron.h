// -*-c++-*-
/// \file perceptron.h
/// \brief interface to CPerceptron
/// \author Artem Sokolov

#ifndef PERCEPTRON_H__INCLUDED
#define PERCEPTRON_H__INCLUDED

#include "types.h"

#include <vector>
#include <memory>

using std::vector;
using std::shared_ptr;

/// A collection of parameters defining the behavior of a CPerceptron object
struct SPerceptronParams
{
  /// Set to true if the negative penalty is the loss; false if it's -1
  bool bLossUpdate;

  /// Desired margin for every fold (values are recycled)
  double margin;
};

/// Kernel perceptron over the structured outputs
template< typename _I, typename _O >
class CPerceptron : public CClassifier< _I, _O >
{
public:
  /// Constructor
  CPerceptron( const SPerceptronParams& pp );

  /// Destructor
  ~CPerceptron() {}

protected:
  /// The associated parameters
  SPerceptronParams params;

  /// The coefficients
  /** Maps (index i into the input space set, index j into the output space set) -> alpha_ij */
  smat_t alpha;

private:
  /// Trains the classifier
  virtual void train();

public:
  /// Clears the learned information for the classifier
  virtual void clear() {alpha.clear();}

  /// Performs an update based on a single example
  /** \return Loss of the prediction
   */
  double singleUpdate( unsigned int dsi, double pen_scale = 1.0 );

public:
  /// Displays the coefficients to an arbitrary output stream
  void display( std::ostream& os = std::cout ) const;

  /// Computes argmax_y f(xi,y) over Y
  /** All indices are over the dataset CPerceptron::pds */
  unsigned int argmax( unsigned int xi,
		       const vector<unsigned int>::const_iterator& Ybegin,
		       const vector<unsigned int>::const_iterator& Yend ) const;

private:
  /// Updates alpha_ij by delta
  void updateAlpha( int i, int j, double delta );

public:
  /// Generic implementation of the compatibility function
  virtual double f( pair< const CDataSet<_I>&, unsigned int > xi, unsigned int dso ) const;

};

#include "perceptron_impl"

#endif
