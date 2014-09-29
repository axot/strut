// go-annotation.cpp - implementation of CGOAnnotation
//
// by Artem Sokolov

#include "go-annotation.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
#include <boost/spirit/include/classic_increment_actor.hpp>
#include <boost/spirit/include/classic_clear_actor.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;

namespace GO
{
    // Returns true if the evidence code is safe
  bool isSafeEvidenceCode( const string& str )
  {
    // Safe evidence codes: 'IDA', 'TAS', 'IMP', 'IGI', 'IPI', 'IEP', 'NAS', 'IC'
    if( str == "IDA" ) return true;
    if( str == "TAS" ) return true;
    if( str == "IMP" ) return true;
    if( str == "IGI" ) return true;
    if( str == "IPI" ) return true;
    if( str == "IEP" ) return true;
    if( str == "NAS" ) return true;
    if( str == "IC" ) return true;
    return false;
  }

  // Constructor
  CGOAnnotation::CGOAnnotation( const string& str )
  {
    using spirit::classic::ch_p;
    using spirit::classic::assign_a;

    if( spirit::classic::parse( str.c_str(),
		       (
			(+(~ch_p('\t')))[assign_a(r_DB)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_ObjectID)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_ObjectSym)] >> ch_p('\t') >>
			/* optional */     (*(~ch_p('\t')))[assign_a(r_Qualifier)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_GOID)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_DBRef)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_ECode)] >> ch_p('\t') >>
			/* optional */     (*(~ch_p('\t')))[assign_a(r_WorF)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_Aspect)] >> ch_p('\t') >>
			/* optional */     (*(~ch_p('\t')))[assign_a(r_ObjectName)] >> ch_p('\t') >>
			/* optional */     (*(~ch_p('\t')))[assign_a(r_Synonym)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_ObjectType)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_Taxon)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_Date)] >> ch_p('\t') >>
			(+(~ch_p('\t')))[assign_a(r_AssignedBy)] >> *(ch_p('\t')) >>
			/* optional */     (*(~ch_p('\t')))[assign_a(r_AnnotExt)] >> *(ch_p('\t')) >>
			/* optional */     (*(~ch_p('\t')))[assign_a(r_GenePFID)] >> *(ch_p('\t'))
			)
		       ).full == false )
      {
	cout<<"Failed to parse "<<str<<endl;
      }
  }

  // Displays the record to the provided output stream
  void CGOAnnotation::display( std::ostream& os ) const
  {
    os<<"DB            : "<<r_DB<<endl;
    os<<"Object ID     : "<<r_ObjectID<<endl;
    os<<"Object Symbol : "<<r_ObjectSym<<endl;
    os<<"Qualifier     : "<<r_Qualifier<<endl;
    os<<"GO ID         : "<<r_GOID<<endl;
    os<<"DB Reference  : "<<r_DBRef<<endl;
    os<<"Evidence Code : "<<r_ECode<<endl;
    os<<"W or F        : "<<r_WorF<<endl;
    os<<"Aspect        : "<<r_Aspect<<endl;
    os<<"Object Name   : "<<r_ObjectName<<endl;
    os<<"Synonym       : "<<r_Synonym<<endl;
    os<<"Object Type   : "<<r_ObjectType<<endl;
    os<<"Taxon         : "<<r_Taxon<<endl;
    os<<"Date          : "<<r_Date<<endl;
    os<<"Assigned By   : "<<r_AssignedBy<<endl;
  }

  // Returns the aspect of the annotation
  ONTOLOGY_INDEX CGOAnnotation::getAspect() const
  {
    if( r_Aspect == "F" ) return GO_MF;
    else if( r_Aspect == "P" ) return GO_BP;
    else if( r_Aspect == "C" ) return GO_CC;
    else return GO_NONE;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Constructor
  CGOACollection::CGOACollection( const string& filename )
  {
    // Open the file, constructing a gzip filter (if necessary)
    boost::iostreams::filtering_istream ifs;
    if( boost::ends_with( filename, ".gz" ) )
      {
	ifs.push(boost::iostreams::gzip_decompressor());
	ifs.push(boost::iostreams::file_source(filename, std::ios_base::binary));
      }
    else
      ifs.push(boost::iostreams::file_source(filename));

    if( ifs.good() == false )
      throw std::runtime_error( "Failed to open a GO annotations file" );

    // Parse the file line-by-line
    string str;
    for( std::getline( ifs, str ); ifs.fail() == false; std::getline( ifs, str ) )
      {
	// Ignore the empty lines
	if( str.length() < 1 ) continue;

	// Ignore the comment lines
	if( str[0] == '!' ) continue;

	// Construct and store the annotation
	CGOAnnotation goa( str );
	annots.push_back( goa );
      }

    postprocessing();
  }

  // Post-processing of raw GO annotation values
  void CGOACollection::postprocessing()
  {
    // For each raw annotation entry...
    for( unsigned int i = 0; i < annots.size(); i++ )
      {
	// Retrieve the object ID
	string objID = annots[i].getObjectID();

	// Search the map to see if an entry with this ID already exists
	if( goa_map.find( objID ) == goa_map.end() )
	  {
	    // First time seeing this ID, compose the entry

	    // Retrieve the object names
	    vector< string > obj_names;
	    obj_names.push_back( annots[i].getObjectSymbol() );

	    // Retrieve the object synonyms
	    string syns = annots[i].getSynonyms();
	    if( syns.length() > 0 )
	      {
		if( spirit::classic::parse( syns.c_str(),
				   (
				    (+(~spirit::classic::ch_p('|')))[spirit::classic::push_back_a(obj_names)] >>
				    *( spirit::classic::ch_p('|') >> 
				       (+(~spirit::classic::ch_p('|')))[spirit::classic::push_back_a(obj_names)] )
				    )
				   ).full == false )
		  throw std::runtime_error( "Failed to process GO annotation names/synonyms" );
	      }

	    // Compose and add the entry
	    vector< unsigned int > i_list; i_list.push_back( i );
	    std::pair< vector< string >, vector< unsigned int > > value( obj_names, i_list );
	    goa_map[objID] = value;
	  }
      
	else
	  {
	    // Already have this ID, simply add the index
	    goa_map[objID].second.push_back(i);
	  }
      }

    // Compose all protein names (including synonyms) into a single structure
    for( goa_map_citer_t iter = goa_map.begin(); iter != goa_map.end(); iter++ )
      {
	// Insert the ObjectID
	prot_names.insert( iter->first );

	// Insert the synonyms
	for( unsigned int i = 0; i < iter->second.first.size(); i++ )
	  prot_names.insert( iter->second.first[i] );
      }
  }

  // Takes a protein name, searches through all synonyms and returns an iterator into goa_map
  // Returns goa_map.end() if no match found
  map< string, std::pair< vector< string >, vector< unsigned int > > >::const_iterator 
  CGOACollection::findProtein( const string& prot_name, bool bCaseSensitive ) const
  {
    using boost::to_upper;
    using boost::to_upper_copy;

    // Attempt to find the name first
    for( goa_map_citer_t iter = goa_map.begin(); iter != goa_map.end(); ++iter )
      {
	if( bCaseSensitive && iter->first == prot_name ) return iter;
	if( !bCaseSensitive && to_upper_copy( iter->first ) == to_upper_copy( prot_name ) ) return iter;
      }

    // If failed, search the synonyms
    for( goa_map_citer_t iter = goa_map.begin(); iter != goa_map.end(); ++iter )
      {
	string s = prot_name;
	if( !bCaseSensitive ) to_upper( s );

	// Attempt to match the Object names
	for( unsigned int i = 0; i < iter->second.first.size(); ++i )
	  {
	    if( bCaseSensitive && s == iter->second.first[i] ) return iter;
	    if( !bCaseSensitive && s == to_upper_copy( iter->second.first[i] ) ) return iter;
	  }
      }
  
    return goa_map.end();
  }

  // Returns true if the j^th annotation is not a NOT, has a safe evidence code and satisfies the filter
  bool CGOACollection::isAnnotationGood( unsigned int j, ONTOLOGY_INDEX filter ) const
  {
    // Check for NOT
    if( annots[j].getQualifier() == string( "NOT" ) ) return false;

    // Check for safe evidence code
    if( isSafeEvidenceCode( annots[j].getEvidenceCode() ) == false ) return false;

    // Check the filter
    if( (annots[j].getAspect() & filter) > 0 )
      return true;
    else
      return false;
  }

  // Displays the collection
  void CGOACollection::display( std::ostream& os ) const
  {
    for( goa_map_citer_t iter = goa_map.begin(); iter != goa_map.end(); iter++ )
      {
	os<<iter->first<<":"<<endl;
	os<<"  ";
	for( unsigned int i = 0; i < iter->second.first.size(); i++ )
	  os<<iter->second.first[i]<<" ";
	os<<endl<<"  ";
	for( unsigned int i = 0; i < iter->second.second.size(); i++ )
	  os<<iter->second.second[i]<<" ";
	os<<endl;
      }
  }

  // Returns the protein ID given its synonym
  const string CGOACollection::getObjectID( const string& syn, bool bCaseSensitive ) const
  {
    goa_map_citer_t iter = findProtein( syn, bCaseSensitive );
    if( iter == goa_map.end() ) return string("");
    else return iter->first;
  }

  // Returns a vector of GO IDs for a protein
  vector< string > CGOACollection::getGOIDs( const string& prot_name,
					     ONTOLOGY_INDEX filter,
					     bool bCaseSensitive ) const
  {
    vector< string > res;

    // Locate the protein in the collection
    goa_map_citer_t iter = findProtein( prot_name, bCaseSensitive );
    if( iter == goa_map.end() ) return res;

    // Retrieve the GO IDs
    for( unsigned int i = 0; i < iter->second.second.size(); i++ )
      {
	unsigned int j = iter->second.second[i];
	if( isAnnotationGood( j, filter ) )
	  res.push_back( annots[j].getGOID() );
      }

    return res;
  }

  // Returns whether or not the protein contains GO IDs that match the filter
  bool CGOACollection::hasGOIDs( const string& prot_name,
				 ONTOLOGY_INDEX filter,
				 bool bCaseSensitive ) const
  {
    // Locate the protein in the collection
    goa_map_citer_t iter = findProtein( prot_name, bCaseSensitive );
    if( iter == goa_map.end() ) return false;

    // Retrieve the GO IDs
    for( unsigned int i = 0; i < iter->second.second.size(); i++ )
      {
	unsigned int j = iter->second.second[i];
	if( isAnnotationGood( j, filter ) )
	  return true;
      }

    return false;
  }

  // Returns all synonym names for all proteins that have a non-NOT annotation that falls into
  //   one of the provided categories
  void CGOACollection::getAnnotatedProteins( ONTOLOGY_INDEX filter,
					     set< string >& res ) const
  {
    // Traverse the (protein, annotation info) entries
    for( goa_map_citer_t iter = goa_map.begin(); iter != goa_map.end(); iter++ )
      {
	if( hasGOIDs( iter->first, filter ) )
	  res.insert( iter->first );
      }
  }

  // Provided for convenience
  vector< string> CGOACollection::getAnnotatedProteins( ONTOLOGY_INDEX filter ) const
  {
    // Retrieve the set of proteins
    set< string > pids;
    getAnnotatedProteins( filter, pids );

    // Convert the set to a vector
    vector< string > res( pids.begin(), pids.end() );
    return res;
  }
}

