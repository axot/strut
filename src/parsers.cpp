/// \file parsers.h
/// \brief implementation of a collection of parsers
/// \author Artem Sokolov

#include "parsers.h"
#include "sample.h"

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
#include <boost/spirit/include/classic_increment_actor.hpp>
#include <boost/spirit/include/classic_clear_actor.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

namespace spirit = boost::spirit::classic;

using spirit::parse;
using spirit::blank_p;
using spirit::ch_p;
using spirit::real_p;
using spirit::assign_a;
using spirit::push_back_a;

shared_ptr< CFeatMap > parseSparseFile( const string& filename,
					CDataSet<CSparseSample>& ds,
					char cSIDSep, char cPairSep, char cFVSep )
{
  // Create a feature map
  shared_ptr< CFeatMap > pfmap( new CFeatMap );
  parseSparseFile( filename, ds, pfmap,
		   cSIDSep, cPairSep, cFVSep );
  return pfmap;
}
  

void parseSparseFile( const string& filename,
		      CDataSet<CSparseSample>& ds,
		      shared_ptr< CFeatMap > pfmap,
		      char cSIDSep, char cPairSep, char cFVSep )
{
  // Open the file
  boost::iostreams::filtering_istream ifs;
  if( openReadFile( filename, ifs ) == false )
    throw std::runtime_error( "Failed to open " + filename );

  // Parse the file line-by-line
  string str;
  for( std::getline( ifs, str ); ifs.fail() == false; std::getline( ifs, str ) )
    {
      // Do nothing in case of comment lines
      if( str[0] == '#' ) continue;

      string sampleID;
      vector<string> keys;
      vector<double> vals;

      // Parse the line
      if( parse ( str.c_str(),
		  ( (+(~ch_p(cSIDSep)))[assign_a(sampleID)] >>
		    !(
		      ch_p(cSIDSep) >>
		      *(
			!ch_p(cPairSep) >>
			(+(~ch_p(cFVSep)))[push_back_a(keys)] >> 
			ch_p(cFVSep) >>
			real_p[push_back_a(vals)] >>
			!ch_p(cPairSep)
			)
		      )
		    )).full == false )
	throw std::logic_error( "Failed to parse "+str );

      // Create a sample
      auto_ptr< CSparseSample > psmpl( new CSparseSample(pfmap) );

      // Convert the features to indices, adding new features as necessary
      for( unsigned int i = 0; i < keys.size(); i++ )
	{
	  // Add the feature as necessary
	  unsigned int fi = pfmap->addFeat( keys[i] );

	  // Add the pair to the sample
	  psmpl->setValue( fi, vals[i] );
	}

      // Add the sample to the dataset
      ds.addSample( sampleID, psmpl );
    }
}


// Parses a file in tab-delimited format
void parseTabDelFile( const string& filename,
		      CDataSet< CSparseSample >& ds )
{
  // Create a feature map
  shared_ptr< CFeatMap > pfmap( new CFeatMap );

  // Open the file
  boost::iostreams::filtering_istream ifs;
  if( openReadFile( filename, ifs ) == false )
    throw std::runtime_error( "Failed to open " + filename );

  // Retrieve the header
  string str;
  std::getline( ifs, str );
  if( ifs.fail() == true )
    throw std::runtime_error( "Failed to read the header" );

  // Parse the header
  vector<string> ids;
  if( parse( str.c_str(),
	     ( 
	      !ch_p('\"') >> *(~ch_p('\t') & ~ch_p('\"')) >> !ch_p('\"') >>
	      *( ch_p('\t') >> 
		 !ch_p('\"') >>
		 (+(~ch_p('\t') & ~ch_p('\"')))[push_back_a(ids)] >>
		 !ch_p('\"') ) >>
	      *( ch_p('\t') )
	      )
	     ).full == false )
    throw std::logic_error( "Failed to parse the header in file "+filename );

  // Add the features to the feature map
  for( unsigned int i = 0; i < ids.size(); i++ )
    pfmap->addFeat( ids[i] );
  if( pfmap->nFeats() != ids.size() )
    throw std::logic_error( "Non-unique feature IDs in file "+filename );

  // Parse the content line-by-line
  for( std::getline( ifs, str ); ifs.fail() == false; std::getline( ifs, str ) )
    {
      string sampleID;
      vector<double> values;

      // Parse the line
      if( parse( str.c_str(),
		 (
		  (+(~ch_p('\t')))[assign_a(sampleID)] >>
		  *( *(blank_p) >> real_p[push_back_a(values)] ) >>
		  *(blank_p)
		  ) ).full == false )
	throw std::logic_error( "Failed to parse "+str );
      if( values.size() != pfmap->nFeats() )
	throw std::logic_error( "Inconsistent file "+str );

      // Create the sample
      auto_ptr< CSparseSample > psmpl( new CSparseSample(pfmap) );

      // Add the values
      for( unsigned int i = 0; i < values.size(); i++ )
	psmpl->setValue( i, values[i] );

      // Add the sample to the dataset
      ds.addSample( sampleID, psmpl );
    }
}

/// Loads a collection of files, one per kernel space
CDataSet<vSparseSample> loadKernels( vector< string >& filenames )
{
  CDataSet<vSparseSample> ds;
  for( unsigned int i = 0; i < filenames.size(); ++i )
    {
      string fn = filenames[i];
      cout<<"Loading "<<fn<<" "; cout.flush();
      CDataSet<CSparseSample> d;
      
      // Infer the format
      if( fn.find( "sdat" ) != fn.npos )
	{
	  cout<<"using .sdat parser ( , , = )   "; cout.flush();
	  parseSparseFile( fn, d, ',', ',', '=' );
	}
      else
	{
	  cout<<"using .dat parser ( ,   : )   "; cout.flush();
	  parseSparseFile( fn, d, ',', ' ', ':' );
	}
      
      cout<<d.size()<<" samples; "<<nFeats( d )<<" features"<<endl;
      
      // Append the new dataset to the multiset
      bool bRemoveMissing = (i > 0);	// Dont remove missing for the first space
      expand( ds, d, bRemoveMissing );
    }
  cout<<"Joint set has "<<ds.size()<<" samples; "<<nKernels( ds )<<" kernels for "<<nFeats( ds )<<" features total"<<endl;
  return ds;
}
