// misc.cpp - miscellaneous
//
// by Artem Sokolov

#include "misc.h"

#include <cmath>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <boost/algorithm/string/predicate.hpp>

// Random number generator
boost::mt19937 g_RNG;

/// Infinity
double getInfinity() { return std::numeric_limits< double >::infinity(); }

// Generic open file
// Determines if the file is gzipped and constructs a filter as appropriate
// Returns true if successful, false otherwise
bool openReadFile( const string& filename, boost::iostreams::filtering_istream& ifs )
{
  // Attempt to open a file to ensure its existence
  std::ifstream file( filename.c_str() );
  if( file.good() == false ) return false;
  file.close();

  // Construct a gzip filter (if necessary)
  if( boost::ends_with( filename, ".gz" ) )
    {
      ifs.push(boost::iostreams::gzip_decompressor());
      ifs.push(boost::iostreams::file_source(filename, std::ios_base::binary));
    }
  else
    ifs.push(boost::iostreams::file_source(filename));

  return true;
}

// Open a gzipped file for writing
void openWriteFile( const string& filename, boost::iostreams::filtering_ostream& ofs )
{
  // Construct a gzip filter
  ofs.push( boost::iostreams::gzip_compressor() );
  ofs.push( boost::iostreams::file_sink(filename, std::ios_base::binary) );
}

