/// \file blastout.cpp
/// \brief implementation of CBLASTOutput
/// \author Artem Sokolov

#include "blastout.h"
#include "sample.h"

#include <sstream>
#include <fstream>

namespace GO
{
  CBLASTOutput::SBOEntry::SBOEntry( std::string raw_line )
  {
    // Use a string stream to parse the tab-delimited entries in the line
    std::istringstream iss( raw_line );

    iss>>query_id;
    iss>>subject_id;
    iss>>percent_identity;
    iss>>alignment_length;
    iss>>mismatches;
    iss>>gaps;
    iss>>q_start;
    iss>>q_end;
    iss>>s_start;
    iss>>s_end;
    iss>>e_value;
    iss>>bit_score;
  }

  void CBLASTOutput::SBOEntry::display( std::ostream& os ) const
  {
    os<<query_id<<"\t";
    os<<subject_id<<"\t";
    os<<percent_identity<<"\t";
    os<<alignment_length<<"\t";
    os<<mismatches<<"\t";
    os<<gaps<<"\t";
    os<<q_start<<"\t";
    os<<q_end<<"\t";
    os<<s_start<<"\t";
    os<<s_end<<"\t";
    os<<e_value<<"\t";
    os<<bit_score<<std::endl;
  }

  // Constructor
  CBLASTOutput::CBLASTOutput( const string& filename )
  {
    // Open the file
    boost::iostreams::filtering_istream ifs;
    if( openReadFile( filename, ifs ) == false )
      throw std::runtime_error( "Failed to open a BLAST output file" );
  
    // Parse it line-by-line
    std::string str;
    for( std::getline( ifs, str ); ifs.fail() == false; std::getline( ifs, str ) )
      {
	SBOEntry entry( str );
	std::string qid = entry.query_id;

	// Check to see if we already have an entry with this query_id
	if( entries.find( qid ) == entries.end() )
	  {
	    // Create a new entry
	    std::vector< SBOEntry > v;
	    v.push_back( entry );
	    entries[qid] = v;
	  }

	else
	  {
	    // Append to the entry
	    entries[qid].push_back( entry );
	  }
      }
  }

  // Displays the blast output to an arbitrary stream
  void CBLASTOutput::display( std::ostream& os ) const
  {
    for( emap_const_iter_t iter = entries.begin(); iter != entries.end(); ++iter )
      for( unsigned int i = 0; i < iter->second.size(); i++ )
	iter->second[i].display( os );
  }

  /// Returns the hits associated with the requested key
  std::vector< CBLASTOutput::SBOEntry > CBLASTOutput::getHits( const std::string& key ) const
  {
    emap_const_iter_t iter = entries.find( key );
    if( iter == entries.end() ) return std::vector< CBLASTOutput::SBOEntry >();
    else return iter->second;
  }

  // Returns the names of all entries
  std::vector< string > CBLASTOutput::getNames() const
  {
    std::vector< string > res;
    for( emap_const_iter_t iter = entries.begin(); iter != entries.end(); ++iter )
      res.push_back( iter->first );
    return res;
  }

  // Returns true if two proteins are within e_thresh of each other
  bool CBLASTOutput::proximityEVal( const string& p1, const string& p2, double e_thresh ) const
  {
    // Search for p2 in the hits of p1
    emap_const_iter_t iter = find( p1 );
    if( iter != end() )
      {
	// Traverse the hits
	for( unsigned int i = 0; i < iter->second.size(); i++ )
	  {
	    if( iter->second[i].subject_id == p2 &&
		iter->second[i].e_value < e_thresh ) return true;
	  }
      }

    // Search for p1 in the hits of p2
    iter = find( p2 );
    if( iter != end() )
      {
	// Traverse the hits
	for( unsigned int i = 0; i < iter->second.size(); i++ )
	  {
	    if( iter->second[i].subject_id == p1 &&
		iter->second[i].e_value < e_thresh ) return true;
	  }
      }

    return false;
  }

  // Returns true if two proteins are within pi_thresh of each other
  bool CBLASTOutput::proximityPIden( const string& p1, const string& p2, double pi_thresh ) const
  {
    // Search for p2 in the hits of p1
    emap_const_iter_t iter = find( p1 );
    if( iter != end() )
      {
	// Traverse the hits
	for( unsigned int i = 0; i < iter->second.size(); i++ )
	  {
	    if( iter->second[i].subject_id == p2 &&
		iter->second[i].percent_identity > pi_thresh ) return true;
	  }
      }

    // Search for p1 in the hits of p2
    iter = find( p2 );
    if( iter != end() )
      {
	// Traverse the hits
	for( unsigned int i = 0; i < iter->second.size(); i++ )
	  {
	    if( iter->second[i].subject_id == p1 &&
		iter->second[i].percent_identity > pi_thresh ) return true;
	  }
      }

    return false;
  }

}
