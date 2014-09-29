// -*-c++-*-
/// \file types.h
/// \brief A collection of typedefs specific to this package
/// \author Artem Sokolov

#ifndef TYPES_H__INCLUDED
#define TYPES_H__INCLUDED

#include <exception>
#include <stdexcept>
#include <unordered_map>
#include <memory>
#include <functional>

#include <string>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <iterator>
#include <limits>

using std::string;
using std::vector;
using std::set;
using std::pair;
using std::map;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::shared_ptr;
using std::function;

/// (int, int) hasher
struct uint_pair_hasher
{
  /// The hash function
  size_t operator() (const std::pair< unsigned int, unsigned int >& key ) const
  {
    return (key.first << ( sizeof(unsigned int) * 8/2)) ^ key.second;
  }
};

/// Displays a vector
template< typename _T >
std::ostream& operator<<( std::ostream& os, const vector< _T >& v )
{ std::copy( v.begin(), v.end(), std::ostream_iterator<_T>( os, " " ) ); return os; }

/// A range of indices
typedef std::pair< unsigned int, unsigned int > irange_t;

/// A vector of ranges
typedef std::vector< irange_t > virange_t;

/// Sparse vector: i -> V_i
typedef std::unordered_map< unsigned int, double > svec_t;

/// Sparse matrix: (i,j) -> M_ij
typedef std::unordered_map< std::pair< unsigned int, unsigned int >, double, uint_pair_hasher > smat_t;

/// A string -> integer map
typedef std::unordered_map< string, int > simap_t;

/// A string -> unsigned integer map
typedef std::unordered_map< string, unsigned int > sumap_t;

/// An unsigned int -> unsigned int map
typedef std::unordered_map< unsigned int, unsigned int > uumap_t;

#endif
