// -*-c++-*-
/// \file clsf.h
/// \brief CGenericClassifier interface
/// \author Artem Sokolov

#ifndef CLSF_H__INCLUDED
#define CLSF_H__INCLUDED

#include "io-dataset.h"
#include "params.h"
#include "infer.h"

#include <memory>

using std::vector;
using std::shared_ptr;

// Forward declaration
template< typename _I, typename _O > class CClassifier;

/// Computes the compatibility function values for xk and a range of labels
template< typename _I, typename _O >
class FFComputer
{
public:
  FFComputer( const CClassifier<_I,_O>& clsf_,
	      unsigned int xk_,
	      const std::pair< unsigned int, unsigned int >& range_,
	      vector< double >& res_ )
    : clsf( clsf_ ), xk( xk_ ), range( range_ ), res( res_ ) {}

private:
  const CClassifier<_I,_O>& clsf;
  unsigned int xk;
  std::pair<unsigned int, unsigned int> range;
  vector< double >& res;

public:
  void operator() ()
  {
    for( unsigned int yi = range.first; yi < range.second; yi++ )
      res[yi] = clsf.f( xk, yi );
  }
};

/// The base classifier class
template< typename _I, typename _O >
class CClassifier
{
private:
  typedef function< unsigned int( pair< const CDataSet<_I>&, unsigned int > ) > inferer_t;
  typedef FLoopyArgmax<_I,_O> defInferer;

public:
  /// Constructor
  explicit CClassifier( string str = string("") )
    : name( str ), inferer( defInferer( *this ) ) {}

  /// Constructor
  explicit CClassifier( string str, const inferer_t& inf )
    : name( str ), inferer( inf ) {}

  /// Destructor
  virtual ~CClassifier() {}

protected:
  /// The associated training data
  shared_ptr< const CIODataSet<_I,_O> > pdsTrain;

  /// The public name of the classifier
  string name;

  /// Used to perform inference
  inferer_t inferer;

public:
  /// Retrieve the public name
  const string getName() const
  { return name; }

  /// Provides a new inference algorithm
  void setInferer( const inferer_t& inf ) { inferer = inf; }

  /// Retrieves the associated dataset
  shared_ptr< const CIODataSet<_I,_O> > getDataSet() const {return pdsTrain;}

  /// Displays a message as pertaining to this classifier
  void displayMessage( string msg, std::ostream& os = std::cout ) const
  { os<<name<<" : "<<msg<<std::endl; }

  /// Trains the classifier on a specific dataset
  void train( shared_ptr< const CIODataSet<_I,_O> > pds );

  /// Makes predictions for a new set of data. Saves predictions to file fnPred
  vector<unsigned int> predict( shared_ptr< const CDataSet<_I> > ds, const string& fnPred = "" ) const;

  /// Tests the classifier on a new set of data
  /** \return Loss values for each sample
      \param[in] dsTest dataset to test the classifier on
      \param[in] fnPred filename for saving predictions to
   */
  vector<double> test( shared_ptr< const CIODataSet<_I, _O> > dsTest, const string& fnPred = "" ) const;

  /// Infers the most compatible label for a foreign sample
  unsigned int infer( const CDataSet<_I>& ds, unsigned int xi ) const
  {
    std::pair< const CDataSet<_I>&, unsigned int > dsx( ds, xi );
    return inferer( dsx );
  }

  /// Infers the most compatible label for a foreign sample
  unsigned int infer( std::pair< const CDataSet<_I>&, unsigned int > dsx ) const
  { return inferer( dsx ); }

  /// Computes the compatibility
  double f( const CDataSet<_I>& ds, unsigned int xi, unsigned int yj ) const
  {
    std::pair< const CDataSet<_I>&, unsigned int > dsx( ds, xi );
    return f( dsx, yj );
  }

  /// Computes the compatibility for an internal sample
  double f( unsigned int xi, unsigned int yj ) const
  {
    std::pair< const CDataSet<_I>&, unsigned int > dsx( *pdsTrain->getI(), xi );
    return f( dsx, yj );
  }

public:
  /// Clear the learned information for the classifier
  virtual void clear() = 0;

  /// Computes the compatibility score for a foreign example
  virtual double f( std::pair< const CDataSet<_I>&, unsigned int > dsx, unsigned int yj ) const = 0;

  /// Saves the classifier to a file
  virtual void save( const string& filename ) const
  { throw std::logic_error( "Function is not used" ); }

  /// Loads the classifier associated with the provided data from a file; returns false if failed
  virtual bool load( shared_ptr< const CIODataSet<_I, _O> > pds, const string& filename )
  { throw std::logic_error( "Function is not used" ); return false; }

private:
  /// Performs inference on a range of test samples
  class CTester
  {
  public:
    CTester( const CClassifier<_I,_O>& clsf_,
	     const CDataSet<_I>& dsTest_,
	     const irange_t& iRange_,
	     const irange_t& oRange_,
	     vector<unsigned int>& res_ )
      : clsf( clsf_ ), dsTest( dsTest_ ), iRange( iRange_ ), oRange( oRange_ ), res( res_ ) {}

  private:
    const CClassifier<_I,_O>& clsf;
    const CDataSet<_I>& dsTest;
    irange_t iRange;
    irange_t oRange;
    vector<unsigned int>& res;

  public:
    void operator() ();
  };

  /// Trains the classifier
  virtual void train() = 0;
};

#include "clsf_impl"

#include "random-clsf.h"
#include "perceptron.h"
#include "nssvm.h"
#include "ssvm.h"

#endif
