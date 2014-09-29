// -*-c++-*-
/// \file misc.h
/// \brief Miscellaneous functions and definitions
/// \author Artem Sokolov

#ifndef MISC_H__INCLUDED
#define MISC_H__INCLUDED

#include "types.h"
#include "kernel.h"

#include <boost/iostreams/filtering_stream.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

/// Random number generator
extern boost::mt19937 g_RNG;

/// Infinity
double getInfinity();

/// Generic open file
/** Determines if the file is gzipped and constructs a filter as appropriate
    \return true if successful, false otherwise
*/
bool openReadFile( const string& filename, boost::iostreams::filtering_istream& ifs );

/// Opens a gzipped file for writing
void openWriteFile( const string& filename, boost::iostreams::filtering_ostream& ofs );

#endif
