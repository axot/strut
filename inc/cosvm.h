// -*-c++-*-
/// \file cosvm.h
/// \brief Interface to CCOSVM
/// \author Artem Sokolov

#ifndef COSVM_H__INCLUDED
#define COSVM_H__INCLUDED

#include <memory>

#include "io-dataset.h"

using std::shared_ptr;

/// The set of parameters associated with CCOSVM
struct SCOSVMParams
{
  /// Penalty for margin violation by the labeled examples (per sample)
  double Cn_l;
  
  /// Penalty for margin violation by the unlabeled examples (per sample)
  double Cn_u;

  /// Desired precision for margin violation
  double eps;

  /// Maximum number of attempts to match up the labels for unlabeled samples
  unsigned int rmax;

  /// true if the unlabeled examples are handled transductively, false if through co-training
  bool bTrans;

  /// Prefix for file saving
  string fnPrefix;
};

/// Two-view support vector machine
template< typename _I1, typename _I2, typename _O >
class CCOSVM
{
private:
  /// The type used to denote the first view
  typedef char _V1;
  
  /// The type used to denote the second view
  typedef int _V2;

  /// A constant of the first-view type
  static const _V1 _v1 = 'f';
  
  /// A constant of the second-view type
  static const _V2 _v2 = 3;

public:
  /// Constructor
  explicit CCOSVM( const SCOSVMParams& p ) :
    params( p ) {}

  /// Destructor
  ~CCOSVM() {}

private:
  /// Trains a particular view on labeled examples
  template< typename _V >
  class CViewTrainer
  {
  public:
    CViewTrainer( _V v_, const string& hdr_, bool bSinglePass_, CCOSVM* pSVM_, 
		  const std::pair< unsigned int, unsigned int >& range_,
		  unsigned int& nAdded_ )
      : v(v_), hdr(hdr_), bSinglePass(bSinglePass_), pSVM(pSVM_), range(range_),
	nAdded( nAdded_ ) {}

  private:
    _V v;
    string hdr;
    bool bSinglePass;
    CCOSVM* pSVM;
    std::pair< unsigned int, unsigned int > range;
    unsigned int& nAdded;
    unsigned int nNew;

  public:
    void operator() ();    
  };

  /// Computes the f values for a range of output indices on internal samples
  template< typename _V >
  class CFComputer
  {
  public:
    CFComputer( _V v_, const CCOSVM* pSVM_, unsigned int xk_,
		const std::pair< unsigned int, unsigned int >& range_,
		vector< double >& res_ )
      : v(v_), pSVM( pSVM_ ), xk( xk_ ), range( range_ ), res( res_ ) {}
    
  private:
    _V v;
    const CCOSVM* pSVM;
    unsigned int xk;
    std::pair< unsigned int, unsigned int > range;
    vector< double >& res;

  public:
    void operator() ();
  };

  /// Computes the f values for a range of output indices on external samples
  template< typename _T1, typename _T2 >
  class CFComputerExt
  {
  public:
    CFComputerExt( const CCOSVM* pSVM_, _T1 x1_, _T2 x2_,
		   const std::pair< unsigned int, unsigned int >& range_,
		   vector< double >& res_ )
      : pSVM( pSVM_ ), x1( x1_), x2( x2_ ), range( range_ ), res( res_ ) {}

  private:
    const CCOSVM* pSVM;
    _T1 x1;
    _T2 x2;
    std::pair< unsigned int, unsigned int > range;
    vector< double >& res;

  public:
    void operator() ()
    {
      for( unsigned int i = range.first; i < range.second; i++ )
	res[i] = pSVM->f( x1, x2, i );
    }
  };

  /// Performs inference on a range of test samples
  class CTester
  {
  public:
    CTester( const CCOSVM* pSVM_,
	     const CIODataSet<_I1, _O>& pTest1_,
	     const CIODataSet<_I2, _O>& pTest2_,
	     const std::pair< unsigned int, unsigned int >& range_,
	     vector<unsigned int>& pred_ )
      : pSVM( pSVM_ ), pTest1( pTest1_ ), pTest2( pTest2_ ), 
	range( range_ ), pred( pred_ ) {}

  private:
    const CCOSVM* pSVM;
    const CIODataSet<_I1, _O>& pTest1;
    const CIODataSet<_I2, _O>& pTest2;
    std::pair< unsigned int, unsigned int > range;
    vector<unsigned int>& pred;

  public:
    void operator() ();
  };

private:
  /// The associated parameters
  SCOSVMParams params;

  /// The associated training data for view 1
  shared_ptr< const CIODataSet<_I1, _O> > pds1;

  /// The associated training data for view 2
  shared_ptr< const CIODataSet<_I2, _O> > pds2;

  /// The view 1 weight vector
  vector< svec_t > alpha1;

  /// The view 2 weight vector
  vector< svec_t > alpha2;

  /// The sum of per-sample weights for view 1
  vector< double > asum1;

  /// The sum of per-sample weights for view 2
  vector< double > asum2;

  /// Assignments to unlabeled examples to be treated as "true" by view 1
  uumap_t yassign1;

  /// Assignments to unlabeled examples to be treated as "true" by view 2
  uumap_t yassign2;

public:
  const shared_ptr< const CDataSet<_O> > getOutputSpace() const
  { return pds1->getO(); }

private:

  /// Adds a new constraint to a particular view; throws if constraint exists
  void addConstraint_impl( vector< svec_t >& alpha, unsigned int xk, unsigned int ybar );

  /// Adds a new constraint to view 1
  void addConstraint( _V1, unsigned int xk, unsigned int ybar )
  { addConstraint_impl( alpha1, xk, ybar ); }

  /// Adds a new constraint to view 2
  void addConstraint( _V2, unsigned int xk, unsigned int ybar )
  { addConstraint_impl( alpha2, xk, ybar ); }

  ///////////////////////////////////////////////////////////////////////

  /// Updates the sum of alphas for a particular sample in a particular view
  void updateSum_impl( const vector< svec_t >& alpha, vector< double >& asum, unsigned int xk );

  /// Updates the sum of a sample in view 1
  void updateSum( _V1, unsigned int xk ) {updateSum_impl( alpha1, asum1, xk ); }

  /// Updates the sum of a sample in view 2
  void updateSum( _V2, unsigned int xk ) {updateSum_impl( alpha2, asum2, xk ); }

  ///////////////////////////////////////////////////////////////////////

  /// Returns the alpha sum for sample xk in view 1
  double getSum( _V1, unsigned int xk ) const {return asum1[xk];}

  /// Returns the alpha sum for sample xk in view 2
  double getSum( _V2, unsigned int xk ) const {return asum2[xk];}

  ///////////////////////////////////////////////////////////////////////
  
  /// Returns the true / assigned label for a sample in a particular view
  template< typename _I >
  unsigned int map_impl( const CIODataSet<_I, _O>& ds,
			 const uumap_t& yassign, unsigned int xk ) const;

  /// Returns the true / assigned label for a sample in view 1
  unsigned int map( _V1, unsigned int xk ) const { return map_impl( *pds1, yassign1, xk ); }

  /// Returns the true / assigned label for a sample in view 1
  unsigned int map( _V2, unsigned int xk ) const { return map_impl( *pds2, yassign2, xk ); }

  ///////////////////////////////////////////////////////////////////////
  
  /// Computes the compatibility function values for a particular view
  template< typename _V, typename _I, typename _T >
  double f_view_impl( _V v, const CIODataSet<_I, _O>& ds,
		      const vector< svec_t >& alpha, _T x, unsigned int y ) const;

  /// Returns the compatibility function value for view 1
  template< typename _T > double f_view( _V1 v, _T x, unsigned int y ) const 
  { return f_view_impl( v, *pds1, alpha1, x, y ); }

  /// Returns the compatibility function value for view 2
  template< typename _T > double f_view( _V2 v, _T x, unsigned int y ) const 
  { return f_view_impl( v, *pds2, alpha2, x, y ); }

  ///////////////////////////////////////////////////////////////////////

  /// Computes J_{i,y,j,ybar} in a particular view
  template< typename _V, typename _I >
  double J_impl( _V v, const CIODataSet<_I, _O>& ds,
		 unsigned int i, unsigned int y,
		 unsigned int j, unsigned int ybar ) const;

  /// Computes J_{i,y,j,ybar} in view 1
  double J( _V1 v, unsigned int i, unsigned int y,
	    unsigned int j, unsigned int ybar ) const
  { return J_impl( v, *pds1, i, y, j, ybar ); }

  /// Computes J_{i,y,j,ybar} in view 2
  double J( _V2 v, unsigned int i, unsigned int y,
	    unsigned int j, unsigned int ybar ) const
  { return J_impl( v, *pds2, i, y, j, ybar ); }

  ///////////////////////////////////////////////////////////////////////

  /// Finds the most violated constrained for sample xk in a particular view
  template< typename _V, typename _I >
  int mostViolConstraint_impl( _V v, const CIODataSet<_I, _O>& ds,
			       const vector< svec_t >& alpha, unsigned int xk ) const;

  /// Finds the most violated constrained for sample xk in view 1
  int mostViolConstraint( _V1 v, unsigned int xk ) const
  { return mostViolConstraint_impl( v, *pds1, alpha1, xk ); }

  /// Finds the most violated constrained for sample xk in view 2
  int mostViolConstraint( _V2 v, unsigned int xk ) const
  { return mostViolConstraint_impl( v, *pds2, alpha2, xk ); }

  ///////////////////////////////////////////////////////////////////////

  /// Performs ascent in the subspace of sample xk in a particular view
  template< typename _V, typename _I >
  void SVMopt_impl( _V v, const CIODataSet<_I, _O>& ds,
		    vector< svec_t >& alpha, unsigned int xk, double Cn );

  /// Performs ascent in the subspace of sample xk in view 1
  void SVMopt( _V1 v, unsigned int xk, double Cn ) 
  { SVMopt_impl( v, *pds1, alpha1, xk, Cn ); }

  /// Performs ascent in the subspace of sample xk in view 2
  void SVMopt( _V2 v, unsigned int xk, double Cn )
  { SVMopt_impl( v, *pds2, alpha2, xk, Cn ); }

  ///////////////////////////////////////////////////////////////////////

  /// Performs labeled optimization on sample xk in a particular view
  /** Returns true if a new constraint has been added
   */
  template< typename _V >
  bool optLab_impl( _V v, vector< svec_t >& alpha, unsigned int xk );

  /// Performs labeled optimization on sample xk in view 1
  bool optLab( _V1 v, unsigned int xk ) { return optLab_impl( v, alpha1, xk ); }

  /// Performs labeled optimization on sample xk in view 1
  bool optLab( _V2 v, unsigned int xk ) { return optLab_impl( v, alpha2, xk ); }

  ///////////////////////////////////////////////////////////////////////

  /// Computes the gradient associated with (xk, ybar) in a particular view
  template< typename _V, typename _I >
  double gradient_impl( _V v, shared_ptr< const CIODataSet<_I, _O> > pds,
			const vector< svec_t >& alpha, unsigned int xk, unsigned int ybar ) const;

  /// Computes the gradient associated with (xk, ybar) in view 1
  double gradient( _V1 v, unsigned int xk, unsigned int ybar ) const
  { return gradient_impl( v, pds1, alpha1, xk, ybar ); }

  /// Computes the gradient associated with (xk, ybar) in view 2
  double gradient( _V2 v, unsigned int xk, unsigned int ybar ) const
  { return gradient_impl( v, pds2, alpha2, xk, ybar ); }

  ///////////////////////////////////////////////////////////////////////

  /// Infers the top label for a sample in a particular view
  template< typename _V, typename _I, typename _T >
  unsigned int infer_impl( _V v, shared_ptr< const CIODataSet<_I, _O> > pds,
			   const vector< svec_t >& alpha, _T x ) const;

  /// Infers the top label for a sample in view 1
  template< typename _T >
  unsigned int infer( _V1 v, _T x ) const
  { return infer_impl( v, pds1, alpha1, x ); }

  /// Infers the top label for a sample in view 2
  template< typename _T >
  unsigned int infer( _V2 v, _T x ) const
  { return infer_impl( v, pds2, alpha2, x ); }

  ///////////////////////////////////////////////////////////////////////

  /// Finds the largest margin violation and adds it to the working set if it's greater than xi
  /** Returns true if a new constraint was added */
  template< typename _V, typename _I >
  bool unlabViol_impl( _V v, shared_ptr< const CIODataSet<_I, _O> > pds,
		       vector< svec_t >& alpha, unsigned int xk );

  /// Finds the largest margin violation in view 1
  bool unlabViol( _V1 v, unsigned int xk )
  { return unlabViol_impl( v, pds1, alpha1, xk ); }

  /// Finds the largest margin violation in view 2
  bool unlabViol( _V2 v, unsigned int xk )
  { return unlabViol_impl( v, pds2, alpha2, xk ); }

  ///////////////////////////////////////////////////////////////////////

  /// Performs unlabeled optimization on sample xk
  bool optUnlab( unsigned int xk );

  /// Performs transductive optimization on sample xk
  bool optUnlabTrans( unsigned int xk );

public:
  /// Saves the classifier to a file
  void save( const string& filename ) const;

  /// Loads the classifier from a file; returns false if failed
  bool load( const string& filename );

  /// Preloads the most recent partial result; returns iteartion index of the loaded result
  int preload();

  /// Computes the joint compatibility function value
  template< typename _T1, typename _T2 >
  double f( _T1 x1, _T2 x2, unsigned int y ) const
  { return f_view( _v1, x1, y ) + f_view( _v2, x2, y ); }

  /// Performs joint inference for a sample
  template< typename _T1, typename _T2 >
  unsigned int infer( _T1 x1, _T2 x2 ) const;

  /// Trains the classifier on a specific set of data
  void train( shared_ptr< const CIODataSet<_I1, _O> > pTrain1,
	      shared_ptr< const CIODataSet<_I2, _O> > pTrain2 );

  /// Tests the classifier on a new set of data
  vector<double> test( shared_ptr< const CIODataSet<_I1, _O> > pTest1,
		       shared_ptr< const CIODataSet<_I2, _O> > pTest2,
		       const string& true_fname,
		       const string& pred_fname ) const;

  /// Returns the mappings inferred by each view for an unlabeled example, given its id
  std::pair< unsigned int, unsigned int > mapUnlab( const string& id ) const;
};

template< typename _I1, typename _I2 >
CDataSet<CSparseSample> predScores( const CCOSVM< _I1, _I2, CSparseSample >& clsf, shared_ptr< const CDataSet<_I1> > pds1, shared_ptr< const CDataSet<_I2> > pds2 );

#include "cosvm_impl"

#endif
