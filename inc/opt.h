// -*-c++-*-
/// \file opt.h
/// \brief Optimization routines
/// \author Artem Sokolov

#ifndef OPT_H__INCLUDED
#define OPT_H__INCLUDED

#include <vector>
#include <stdexcept>

using std::pair;

/// Solves LDL^T x = b
template< typename _T >
vector<_T> backSubst( const vector<_T>& L,
		      const vector<_T>& LD,
		      const vector<_T>& b )
{
  const unsigned int n = b.size();
  
  // Perform back-substitution to solve LD(y) = b
  vector< _T > y( n );
  y[0] = b[0] / LD[0];
  for( unsigned int i = 1; i < n; i++ )
    {
      _T prev_sum = 0.0;
      const unsigned int r = i*n;
      for( unsigned int j = 0; j < i; j++ )
	prev_sum += LD[r+j] * y[j];
      y[i] = (b[i] - prev_sum) / LD[r+i];
    }

  // Solve (L^T) x = y
  vector< _T > x( n );
  x[n-1] = y[n-1];		       	// L has 1s on the diagonal
  for( int i = n-2; i >= 0; --i )
    {
      _T prev_sum = 0.0;
      for( unsigned int j = i+1; j < n; ++j )
	prev_sum += L[j*n+i] * x[j];	// Using transpose of L
      x[i] = y[i] - prev_sum;		// L has 1s on the diagonal
    }

  return x;
}

// Solves Qx = b and Qy = 1 for a positive semi-definite matrix Q
template< typename _T >
void solveCholesky( const vector<_T>& Q,
		    const vector<_T>& b,
		    vector<_T>& x,
		    vector<_T>& y )
{
  const unsigned int n = b.size();

  // Perform Cholesky decomposition
  vector< _T > L( n*n );
  vector< _T > D( n );

  L[0] = 1.0; D[0] = Q[0];
  for( unsigned int i = 1; i < n; i++ )
    {
      const unsigned int r = i*n;

      // Compute the entries of L
      for( unsigned int j = 0; j < i; j++ )
	{
	  const unsigned int ij = r+j;
	  const unsigned int c = j*n;
	  L[ij] = Q[ij];
	  for( unsigned int k = 0; k < j; k++ )
	    L[ij] -= L[r+k] * L[c+k] * D[k];
	  L[ij] /= D[j];
	}
      
      const unsigned int ii = r+i;
      L[ii] = 1.0;

      // Compute the entry of D
      D[i] = Q[ii];
      for( unsigned int k = 0; k < i; k++ )
	{
	  const double d = L[r+k];
	  D[i] -= d * d * D[k];
	}

      // Handle singular case through recursion
      if( D[i] == 0 )
	{
	  // Remove the i^th row and column
	  vector< _T > Qr( (n-1)*(n-1) );
	  vector< _T > br( n-1 );
	  unsigned int Qr_i = 0;
	  unsigned int br_i = 0;
	  for( unsigned int ri = 0; ri < n; ri++ )
	    {
	      if( ri == i ) continue;

	      // Copy the ri^th row of Q
	      const unsigned int rr = ri*n;
	      for( unsigned int rj = 0; rj < n; rj++ )
		{
		  if( rj == i ) continue;
		  Qr[Qr_i++] = Q[rr+rj];
		}

	      // Copy the ri^th entry of br
	      br[br_i++] = b[ri];
	    }

	  // Recurse
	  solveCholesky( Qr, br, x, y );

	  // Insert 0s in the appropriate places
	  typename vector<_T>::iterator xbeg = x.begin();
	  x.insert( xbeg+i, 0.0 );
	  typename vector<_T>::iterator ybeg = y.begin();
	  y.insert( ybeg+i, 0.0 );

	  return;
	}
    }

  // Compute the product LD
  vector< _T > LD( n*n );
  for( unsigned int i = 0; i < n; i++ )
    {
      const unsigned int r = i*n;
      for( unsigned int j = 0; j <= i; j++ )
	LD[r+j] = L[r+j] * D[j];
    }

  // Solve LDL^T x = b
  x = backSubst( L, LD, b );

  // Solve LDL^T y = (b+1)
  vector< _T > ones( n, 1.0 );
  //  vector< _T > b1( b.begin(), b.end() );
  //  for( unsigned int i = 0; i < n; ++i ) b1[i] = b1[i] + 1;
  y = backSubst( L, LD, ones );
}

/// Maximizes -0.5 x^T Q x + b^T x subject to \sum x <= Cn
/** Returns the best point on \sum x = Cn if any of coefficients are negative
 */
template< typename _T >
vector<_T> weakQuadraticOpt( const vector<_T>& Q,
			     const vector<_T>& b,
			     const _T& Cn )
{
  // Handle one-variable cases
  if( Q.size() == 1 )
    {

      // Handle coefficient of 0
      if( Q[0] == 0 )
	{
	  if( b[0] == 0 ) return vector<_T>( 1, 0 );	// For maximum flexibility of the other coefficients
	  else throw std::logic_error( "System cannot be solved" );
	}

      // Compute the coefficient
      _T val = b[0] / Q[0];
      if( val > Cn ) val = Cn;
      if( val < 0 ) val = 0;
      return vector<_T>( 1, val );
    }

  const unsigned int n = b.size();

  // Determine the unconstrained solution and the descent vector
  vector<_T> x;
  vector<_T> y;
  solveCholesky( Q, b, x, y );

  // Check if there are negative entries
  bool bHasNeg = false;
  for( unsigned int i = 0; i < n; i++ )
    if( x[i] < 0.0 ) bHasNeg = true;

  // Compute the solution sum
  _T xsum = 0.0;
  for( unsigned int i = 0; i < n; i++ ) xsum += x[i];

  // Check to see if we're within bounds
  if( bHasNeg == false && xsum <= Cn )
    return x;

  // Determine the amount of descent
  _T ysum = 0.0;
  for( unsigned int i = 0; i < n; i++ ) ysum += y[i];

  _T lambda = (xsum - Cn);
  if( lambda < 0.0 ) lambda = 0.0;
  lambda /= ysum;

  // Perform the descent
  for( unsigned int i = 0; i < n; i++ )
    x[i] -= lambda * y[i];

  return x;
}

/// Maximizes -0.5 x^T Q x + b^T x subject to \sum x <= Cn and x_i >= 0
template< typename _T >
vector<_T> strongQuadraticOpt( const vector<_T>& Q,
			       const vector<_T>& b,
			       const _T& Cn )
{
  const unsigned int n = b.size();

  // Find the weak optimization solution
  vector< _T > x = weakQuadraticOpt( Q, b, Cn );

  // Find the smallest entry
  unsigned int argmin = 0; _T min = x[0];
  for( unsigned int i = 1; i < n; i++ )
    if( x[i] < min )
      {
	min = x[i];
	argmin = i;
      }

  // Determine if the solution already satisfies the constraints
  if( min >= 0 ) return x;

  // Eliminate the most negative entry by composing the new Q and b
  unsigned int Qi = 0, bi = 0;
  vector<_T> Q2( (n-1)*(n-1) );
  vector<_T> b2( n-1 );
  for( unsigned int i = 0; i < n; i++ )
    {
      // Skip the row corresponding to the smallest entry
      if( i == argmin ) continue;

      // Copy the entry of b
      b2[bi++] = b[i];

      // Copy the row of Q
      const unsigned int r = i*n;
      for( unsigned int j = 0; j < n; j++ )
	{
	  if( j == argmin ) continue;
	  Q2[Qi++] = Q[r+j];
	}
    }

  // Recurse
  x = strongQuadraticOpt( Q2, b2, Cn );
  
  // Reinstall the removed coefficient as 0
  typename vector<_T>::iterator xbeg = x.begin();
  x.insert( xbeg+argmin, 0 );

  return x;
}

/// Selects a pair of variables for the next SMO iteration
template< typename _T >
pair< int, int > pickVars( const vector<_T>& Q,
			   const vector<_T>& b,
			   const vector<_T>& x,
			   const _T& C )
{
  // Determine the dimensionality
  const unsigned int n = b.size();

  // Compute the gradient: -Qx + b
  vector<_T> g( n );
  for( unsigned int i = 0; i < n; i++ )
    {
      g[i] = 0.0;
      const unsigned int r = i*n;

      // -Qx
      for( unsigned int k = 0; k < n; k++ )
	g[i] -= Q[r+k] * x[k];
      
      // + b
      g[i] += b[i];
    }

  //  cout<<"g: ";
  //  for( unsigned int i = 0; i < n; i++ )
  //    cout<<g[i]<<" ";
  //  cout<<endl;

  // Find the maximal violating pair
  int xmax = -1; _T gmax = -1.0 * std::numeric_limits<double>::infinity();
  int xmin = -1; _T gmin =  1.0 * std::numeric_limits<double>::infinity();

  for( unsigned int i = 0; i < n; i++ )
    {
      // Compare to the largest derivative
      if( g[i] > gmax )
	{
	  if( (g[i] > 0.0 && x[i] < C ) ||
	      (g[i] < 0.0 && x[i] > 0 ) )
	    { xmax = i; gmax = g[i]; }
	}

      // Compare to the smallest derivative
      if( g[i] <= gmin )
	{
	  if( (g[i] > 0.0 && x[i] < C ) ||
	      (g[i] < 0.0 && x[i] > 0 ) )
	    { xmin = i; gmin = g[i]; }
	}
    }

  // Handle the same pick twice
  if( xmin == xmax ) xmax = -1;

  return pair< int, int >( xmin, xmax );
}

/// Maximizes -0.5 x^T Q x + b^T x subject to 0 <= x_i <= C
template< typename _T >
vector<_T> quadOptBox( const vector<_T>& Q,
		       const vector<_T>& b,
		       const _T& C,
		       const _T& epsilon,
		       unsigned int nSteps )
{
  // Determine the dimensionality
  const unsigned int n = b.size();

  // Verify Q dimensionality
  if( Q.size() != n*n ) throw std::logic_error( "quadOptBox: Q and b dimensionality mismatching" );

  // Ensure that diagonal entries are non-0
  for( unsigned int i = 0; i < n; i++ )
    {
      if( Q[i*n+i] == 0 )
	throw std::logic_error( "quadOptBox: Q diagonal entires must not be 0" );
    }

  // Ensure the matrix is symmetric
  for( unsigned int i = 0; i < n; i++ )
    for( unsigned int j = 0; j < n; j++ )
      {
	if( Q[i*n+j] != Q[j*n+i] )
	  throw std::logic_error( "quadOptBox: Q must be symmetric" );
      }

  //  cout<<"Maximizing -0.5 x^T Q x + b^T x, subject to 0 <= x_i <= "<<C<<endl;
  /*  cout<<"Q: ";
  for( unsigned int i = 0; i < n*n; i++ ) cout<<Q[i]<<" ";
  cout<<endl<<"b: ";
  for( unsigned int i = 0; i < n; i++ ) cout<<b[i]<<" ";
  cout<<endl;*/

  // Start with a corner solution
  vector<_T> x( n, 0.0 );

  //  while( true )
  for( unsigned int iter = 0; iter < nSteps; iter++ )
    {
      //      cout<<"x:";
      //      for( unsigned int i = 0; i < n; i++ )
      //	cout<<" "<<x[i];
      //      cout<<endl;
      
      //      cout<<endl;
      //      cout<<"Picking variables"<<endl;
      // Find the maximal violating pair
      pair< int, int > violPair = pickVars( Q, b, x, C );
      //      cout<<"Picked "<<violPair.first<<" "<<violPair.second<<endl;

      // Check the terminating criterion
      if( violPair.first < 0 && violPair.second < 0 ) break;

      // Single-direction optimization
      else if( violPair.first < 0 || violPair.second < 0 )
	{
	  // Retrieve the dimension index
	  unsigned int i;
	  if( violPair.first < 0 ) i = violPair.second;
	  if( violPair.second < 0 ) i = violPair.first;

	  // Save the current value
	  _T oldVal = x[i];

	  // Compute the numerator
	  _T num = b[i];
	  for( unsigned int k = 0; k < n; k++ )
	    num -= Q[i*n+k] * x[k];
	  num += Q[i*n+i] * x[i];

	  // Perform the optimization
	  x[i] = num / Q[i*n+i];

	  // Check for boundaries
	  if( x[i] < 0.0 ) x[i] = 0.0;
	  if( x[i] > C ) x[i] = C;

	  // Compare to the old value and terminate if necessary
	  _T diff = x[i] - oldVal;
	  if( diff < epsilon && diff > -epsilon ) break;
	}

      // Cutting-plane optimization
      else
	{
	  // Fetch the dimension indices
	  unsigned int i = violPair.first;
	  unsigned int j = violPair.second;

	  // Fetch the required entries of Q
	  const _T& qii = Q[i*n+i];
	  const _T& qij = Q[i*n+j];
	  const _T& qjj = Q[j*n+j];

	  // Save the current values
	  _T oldVali = x[i];
	  _T oldValj = x[j];

	  // Compute the right-hand side
	  _T rhsi = b[i]; _T rhsj = b[j];
	  for( unsigned int k = 0; k < n; k++ )
	    { rhsi -= Q[i*n+k] * x[k];
	      rhsj -= Q[j*n+k] * x[k]; }
	  rhsi += qii * x[i] + qij * x[j];
	  rhsj += qij * x[i] + qjj * x[j];

	  // The dimensions are not correlated; deal with each independently
	  if( qij == 0 )
	    {
	      x[i] = rhsi / qii;
	      x[j] = rhsj / qjj;
	    }
	  
	  // Solve a linear system with two equations and two variables
	  else
	    {
	      x[j] = ( qij*rhsi - qii*rhsj ) / ( qij*qij - qii*qjj );
	      x[i] = ( rhsi - qij*x[j] ) / qii;
	    }

	  // Check for boundaries
	  if( x[i] < 0.0 ) x[i] = 0.0; if( x[i] > C ) x[i] = C;
	  if( x[j] < 0.0 ) x[j] = 0.0; if( x[j] > C ) x[j] = C;

	  // Compare to the old value and terminate if necessary
	  _T diffi = x[i] - oldVali;
	  _T diffj = x[j] - oldValj;
	  if( diffi < epsilon && diffi > -epsilon && diffj < epsilon && diffj > -epsilon )
	    break;
	}
    }

  /*
  cout<<"x:";
  for( unsigned int i = 0; i < x.size(); i++ )
    cout<<" "<<x[i];
  cout<<endl;

  // Compute the gradient: -Qx + b
  vector<_T> g( n );
  for( unsigned int i = 0; i < n; i++ )
    {
      g[i] = 0.0;
      const unsigned int r = i*n;

      // -Qx
      for( unsigned int k = 0; k < n; k++ )
	g[i] -= Q[r+k] * x[k];
      
      // + b
      g[i] += b[i];
    }

  cout<<"g:";
  for( unsigned int i = 0; i < g.size(); i++ )
    cout<<" "<<g[i];
  cout<<endl;
  */

  return x;
}

#endif
