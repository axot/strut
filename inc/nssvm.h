// -*-c++-*-
/// \file nssvm.h
/// \brief Interface to CnsSSVM
/// \author Artem Sokolov

#ifndef NSSVM_H__INCLUDED
#define NSSVM_H__INCLUDED

#include <memory>
#include "opt.h"

using std::shared_ptr;

template< typename _I, typename _O, char _R > class CnsSSVM;

/// A collection of parameters defining the behavior of a CnsSSVM object
struct SSSVMParams
{
  /// Slack penalty coefficient for every fold (values are recycled)
  double Cn;

  /// Desired precision
  double eps;

  /// Maximum number of steps QP solver is allowed to take
  unsigned int nMaxQPSteps;

  /// Prefix for file saving
  string fnPrefix;
};

/// A margin-rescaling oracle that loops over all examples in the output space
template< typename _I, typename _O >
struct FLoopyOracle_m
{
  /// The associated classifier
  const CnsSSVM<_I, _O, 'm'>& svm;

  /// The desired precision
  double epsilon;

  /// Constructor
  FLoopyOracle_m( const CnsSSVM<_I, _O, 'm'>& svm_, double eps ) 
    : svm( svm_ ), epsilon( eps ) {}

  /// Returns the most violated constraint for example xi, or -1 if there are none
  int operator() (unsigned int xi);
};

/// A slack-rescaling oracle that loops over all examples in the output space
template< typename _I, typename _O >
struct FLoopyOracle_s
{
  /// The associated classifier
  const CnsSSVM<_I, _O, 's'>& svm;

  /// The desired precision
  double epsilon;

  /// Constructor
  FLoopyOracle_s( const CnsSSVM<_I, _O, 's'>& svm_, double eps ) 
    : svm( svm_ ), epsilon( eps ) {}

  /// Returns the most violated constraint for example xi, or -1 if there are none
  int operator() (unsigned int xi) {throw std::logic_error( "Not yet implemented" );}
};

/// n-slack formulation of the Structured Support Vector Machine
template< typename _I, typename _O, char _R >
class CnsSSVM : public CClassifier<_I, _O>
{
private:
  /// Oracle type (takes an index into the training set, returns index to the most violated constraint, or -1 if there are none)
  typedef function< int( unsigned int ) > oracle_t;

  /// Used to distinguish among template arguments
  template< char _RR > struct _SResc {};

public:
  /// Constructor
  explicit CnsSSVM( const SSSVMParams& pp );

  /// Destructor
  virtual ~CnsSSVM() {}

  /// Creates a default oracle for margin rescaling
  oracle_t makeDefaultOracle( _SResc<'m'> )
  { oracle_t orcl = FLoopyOracle_m<_I,_O>( *this, params.eps ); return orcl; }

  /// Creates a default oracle for slack rescaling
  oracle_t makeDefaultOracle( _SResc<'s'> )
  { oracle_t orcl = FLoopyOracle_s<_I,_O>( *this, params.eps ); return orcl; }

protected:
  /// The associated parameters
  SSSVMParams params;

  /// The separation oracle
  oracle_t oracle;

  /// The coefficients
  vector< svec_t > alpha;

  /// The sum of per-sample coefficients
  vector< double > asum;

  /// Chronological order of the added indices
  vector< vector< unsigned int > > yLast;

  /// Cache of J matrices
  vector< vector< double > > Jcache;

public:
  /// Returns the support vectors associated with example xk
  const svec_t getAlpha( unsigned int xk ) const { return alpha[xk]; }

  /// Trains the classifier
  virtual void train();

  /// Clear the learned information for the classifier
  virtual void clear() {alpha.clear(); asum.clear(); yLast.clear(); Jcache.clear();}

  /// Generic implementation of the compatibility function
  virtual double f( pair< const CDataSet<_I>&, unsigned int > x, unsigned int yj ) const;

  /// Saves the classifier to a file
  virtual void save( const string& filename ) const;

  /// Loads the classifier associated with the provided data from a file; returns false if failed
  virtual bool load( shared_ptr< const CIODataSet<_I, _O> > pds, const string& filename );

  /// Preloads the most recent partial result; returns iteartion index of the loaded result
  int preload( shared_ptr< const CIODataSet<_I, _O> > pds );

private:
  /// Computes the entry of matrix J associated with i, j, y and ybar
  double J( unsigned int i, unsigned int y, unsigned int j, unsigned int ybar ) const;

  /// Computes the dot product between the {k,i}^th row of J and the alpha vector
  double JdotAlpha( unsigned int k, unsigned int i ) const;

  /// Computes the quadratic form of alpha in J: alpha^T J alpha
  double alphaJAlpha() const;

  /// Updates the sum of alphas for the training examples k
  void updateAlphaSum( unsigned int k );

private:
  /// Solves the dual optimization problem using the working set alpha[k]
  void SVMopt( _SResc<'m'>, unsigned int k, double Cn );

  /// Computes the value of the primal and dual objective functions using the current value of the alpha vector
  double obj( _SResc<'m'>, double& dual ) const;

  /// Computes the slack associated with alpha coefficients for the k^th training sample
  /** Uses all constraints / output labels
   */
  double compSlackAll( _SResc<'m'>, unsigned int xi ) const;

  /// Solves the dual optimization problem using the working set alpha[k]
  void SVMopt( _SResc<'s'>, unsigned int k, double Cn )
  { throw std::logic_error( "Not yet implemented" ); }

  /// Computes the value of the primal and dual objective functions using the current value of the alpha vector
  double obj( _SResc<'s'>, double& dual ) const
  { throw std::logic_error( "Not yet implemented" ); }

  /// Finds a label associated with the most violated constraint for sample i
  /** Returns -1 if no constraint is violated by more than params.eps
   */
  int mostViolConstraint( _SResc<'s'>, unsigned int xi ) const
  { throw std::logic_error( "Not yet implemented" ); }

  /// Computes the slack associated with alpha coefficients for the k^th training sample
  /** Uses all constraints / output labels
   */
  double compSlackAll( _SResc<'s'>, unsigned int xi ) const
  { throw std::logic_error( "Not yet implemented" ); }
};

#include "nssvm_impl"

#endif
