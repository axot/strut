// -*-c++-*-
/// \file ssvm.h
/// \brief Interface to C1sSSVM
/// \author Artem Sokolov

#ifndef SSVM_H__INCLUDED
#define SSVM_H__INCLUDED

#include <memory>

using std::shared_ptr;

/// 1-slack formulation of the Structured Support Vector Machine
template< typename _I, typename _O >
class C1sSSVM : public CClassifier<_I, _O >
{ 
public:
  /// Constructor
  explicit C1sSSVM( const SSSVMParams& pp );

  /// Destructor
  virtual ~C1sSSVM() {}

protected:
  /// The associated parameters
  SSSVMParams params;

  /// An assignment of labels to every training sample
  typedef vector<unsigned int> metalabel;

  /// The coefficients
  vector<double> alpha;

  /// The meta-labels comprising the working set
  vector<metalabel> ybar;

  /// Cache of the H matrix
  vector< vector<double> > Hmat;

private:
  /// Computes the most-violated constraints for samples xi_start through xi_end
  class CViolSampleRange
  {
  public:
    /// Constructor
    CViolSampleRange( const C1sSSVM& svm_,
		      unsigned int xk_start, unsigned int xk_end,
		      vector<unsigned int>& res_save ) 
      : svm(svm_), xks(xk_start), xke(xk_end), res(res_save) {}

    /// Destructor
    ~CViolSampleRange() {}

  private:
    /// Pointer to the associated C1sSSVM object
    const C1sSSVM& svm;

    /// Most-violated constraints are computed starting at this sample index
    unsigned int xks;

    /// Most-violated constraints are computed ending at this sample index
    unsigned int xke;

    /// Iterator pointing to where the results should be stored
    vector< unsigned int >& res;

  public:
    /// Executes the search for most-violated constraints
    void operator() ();
  };

public:
  /// Clear the learned information for the classifier
  virtual void clear() {alpha.clear(); ybar.clear(); Hmat.clear();}

  /// Generic implementation of the compatibility function
  virtual double f( pair< const CDataSet<_I>&, unsigned int > x, unsigned int yj ) const;
  
private:
  /// Trains the classifier  
  virtual void train();

  /// Computes the value of xi using the current working set
  double violWS() const;

  /// Computes the value of xi using all constraints / labels
  double violAll() const;

  /// Computes the value of xi associated with a metalabel
  double xi( metalabel lbl ) const;

  /// Computes the partial value of xi associated with a (x_i, ybar)
  double xi( unsigned int x_i, unsigned int y_j ) const;

  /// Generates a random metalabel and adds it to C1sSSVM::ybar
  metalabel makeRndMetalabel();

  /// Generates a metalabel associated with the most violated constraints
  /** Adds the metalabel to C1sSSVM::ybar
      \param[out] amt_tot the total amount of margin violation by the meta-label
   */
  metalabel makeViolMetalabel();

  /// Solves the dual optimization problem
  void SVMopt();

  /// Computes the i^th entry of the gradient
  double del( unsigned int i ) const
  { return (HdotAlpha( i ) - Delta( i )); }

  /// Computes the i^th value of Delta
  /** \param[in] i index into C1sSSVM::ybar */
  double Delta( unsigned int i ) const;

  /// Computes the dot product between the k^th row of H and the alpha vector
  double HdotAlpha( unsigned int k ) const;

  /// Computes and displays the value of the primal and dual objective functions using the current value of the alpha vector
  void displayObj() const;

  /// Caches the H matrix
  void cacheH();

  /// Returns a contribution to the compatibility value based on a kernel and loss values
  virtual double f_contr( double K1, double K2, double LY ) const = 0;

  /// Computes the (i,j)^th value of matrix H
  /** \param[in] i index into C1sSSVM::ybar
      \param[in] j index into C1sSSVM::ybar
  */
  virtual double H( unsigned int i, unsigned int j ) const = 0;

  /// Finds a label associated with the most violated constraint for sample i
  /**
     \param[in] xi an index into the input space dataset
     \param[out] amt the amount of violation
     \return an index into the output space dataset
   */
  virtual unsigned int computeYHat( unsigned int xi ) const = 0;
};

/// An SSVM variant with margin rescaling
template< typename _I, typename _O >
class C1sSSVM_margin : public C1sSSVM<_I, _O>
{
public:
  /// Constructor
  explicit C1sSSVM_margin( const SSSVMParams& pp )
    : C1sSSVM<_I,_O>( pp )
  { this->name = string("1SSVMm"); }

private:
  /// Returns a contribution to the compatibility value based on a kernel and loss values
  virtual double f_contr( double K1, double K2, double LY ) const
  { return (K1 - K2);}

  /// Computes the (i,j)^th value of matrix H
  /** \param[in] i index into C1sSSVM::ybar
      \param[in] j index into C1sSSVM::ybar
  */
  virtual double H( unsigned int i, unsigned int j ) const;

  /// Finds a label associated with the most violated constraint for sample i
  /**
     \param[in] xi an index into the input space dataset
     \param[out] amt the amount of violation
     \return an index into the output space dataset
   */
  virtual unsigned int computeYHat( unsigned int xi ) const;
  
};

/// An SSVM variant with slack rescaling
template< typename _I, typename _O >
class C1sSSVM_slack : public C1sSSVM<_I, _O>
{
public:
  /// Constructor
  explicit C1sSSVM_slack( const SSSVMParams& pp )
    : C1sSSVM<_I,_O>( pp )
  { this->name = "SSVMs"; }

private:
  /// Returns a contribution to the compatibility value based on a kernel and loss values
  virtual double f_contr( double K1, double K2, double LY ) const
  { return LY * (K1 - K2);}

  /// Computes the (i,j)^th value of matrix H
  /** \param[in] i index into C1sSSVM::ybar
      \param[in] j index into C1sSSVM::ybar
  */
  virtual double H( unsigned int i, unsigned int j ) const;

  /// Finds a label associated with the most violated constraint for sample i
  /**
     \param[in] xi an index into the input space dataset
     \param[out] amt the amount of violation
     \return an index into the output space dataset
   */
  virtual unsigned int computeYHat( unsigned int xi ) const;
  
};

#include "ssvm_impl"

#endif
