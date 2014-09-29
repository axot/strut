// go-container.cpp - implementation of CGOContainer
//
// by Artem Sokolov

#include "go-container.h"

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
#include <boost/spirit/include/classic_increment_actor.hpp>
#include <boost/spirit/include/classic_clear_actor.hpp>

using namespace boost;

namespace GO
{
  const ONTOLOGY_INDEX GO_NONE = 0x00;
  const ONTOLOGY_INDEX GO_MF = 0x01;		// Molecular function
  const ONTOLOGY_INDEX GO_BP = 0x02;		// Biological process
  const ONTOLOGY_INDEX GO_CC = 0x04;		// Cellular component

  bool hasMF( ONTOLOGY_INDEX i ) { return (i & GO_MF) == GO_MF; }
  bool hasBP( ONTOLOGY_INDEX i ) { return (i & GO_BP) == GO_BP; }
  bool hasCC( ONTOLOGY_INDEX i ) { return (i & GO_CC) == GO_CC; }

  // Some static definitions
  vector<string>& header_tags()
  {
    static vector<string> header_tags_;
    return header_tags_;
  }

  vector<string>& header_vals()
  {
    static vector<string> header_vals_;
    return header_vals_;
  }

  vector<string>& stanza_names()
  {
    static vector<string> stanza_names_;
    return stanza_names_;
  }

  vector< vector< string > >& stanza_tags()
  {
    static vector< vector< string > > stanza_tags_;
    return stanza_tags_;
  }

  vector< vector< string > >& stanza_vals()
  {
    static vector< vector< string > > stanza_vals_;
    return stanza_vals_;
  }

  // The OBO file format grammar
  class COboGrammar : public spirit::classic::grammar<COboGrammar>
  {
  public:
    template <typename ScannerT>
    class definition
    {
    public:
      definition(const COboGrammar& self)
      {
	// An OBO file consists of a header followed by any number of blank lines, followed by any number stanza
	// The stanzas are also separated by any number of blank lines
	file = header >> *(spirit::classic::ch_p('\n'))
		      >> *(stanza >> *(spirit::classic::ch_p('\n')));

	// The header consists of any number of tag-value pairs
	header = *(tag_value_pair_h >> spirit::classic::ch_p('\n'));

	// A stanza consists of a stanza name followed by any number of tag-value pairs
	stanza = (stanza_name[spirit::classic::push_back_a(stanza_names())] 
		  >> spirit::classic::ch_p('\n') >> *(tag_value_pair_s >> spirit::classic::ch_p('\n')))
	  [spirit::classic::push_back_a(stanza_tags(), cur_stanza_tags)]
	  [spirit::classic::push_back_a(stanza_vals(), cur_stanza_vals)]
	  [spirit::classic::clear_a( cur_stanza_tags )]
	  [spirit::classic::clear_a( cur_stanza_vals )];

	// Stanza name is one of [Term], [Typedef], [Instance]
	stanza_name = spirit::classic::str_p("[Term]") | spirit::classic::str_p("[Typedef]") | spirit::classic::str_p("[Instance]");

	// A tag-value pair is a tag followed by a colon, followed by a value
	tag_value_pair_h = tag[spirit::classic::push_back_a(header_tags())]
	  >> spirit::classic::str_p(": ") >> value[spirit::classic::push_back_a(header_vals())];
	tag_value_pair_s = tag[spirit::classic::push_back_a(cur_stanza_tags)] 
	  >> spirit::classic::str_p(": ") >> value[spirit::classic::push_back_a(cur_stanza_vals)];

	// A tag is a sequence of alphanumeric characters, - and _
	tag = +(spirit::classic::alnum_p | spirit::classic::ch_p('_') | spirit::classic::ch_p('-'));

	// A value is a sequence of any characters
	value = +(~spirit::classic::ch_p('\n'));
      }
    
      spirit::classic::rule<ScannerT> file, header, stanza, stanza_name;
      spirit::classic::rule<ScannerT> tag_value_pair_h, tag_value_pair_s;
      spirit::classic::rule<ScannerT> tag, value;

      spirit::classic::rule<ScannerT> const&
      start() const { return file; }

      vector< string > cur_stanza_tags;
      vector< string > cur_stanza_vals;
    };

    // Constructor
  public:
    COboGrammar() {}
  };

  // Parses an .obo file
  // Returns false in case of a failure, true otherwise
  // Set pof to true if part_of relationships should be counted as is_a for edge construction
  bool CGOContainer::parseOBO( const char* filename, bool pof )
  {
    // Create the file iterators
    spirit::classic::file_iterator<> first( filename );

    if (!first)
      throw std::runtime_error( "Unabled to open " + string(filename) );

    spirit::classic::file_iterator<> last = first.make_end();

    // Parse the file
    COboGrammar g;
    cout<<"Parsing "<<filename<<"... "; cout.flush();
    if( spirit::classic::parse<>( first, last, g ).full )
      cout<<"Success"<<endl;
    else
      {cout<<"Failure"<<endl; return false;}

    // Generate the ontologies
    cout<<"Generating the ontologies"<<endl;
 
    // Temporarily stores directional edges as GOID pairs
    std::multimap< string, string > edges;

    // Traverse the stanzas and generate vertices
    for( unsigned int i = 0; i < stanza_names().size(); i++ )
      {
	// We're only concerned with [Term]s
	if( stanza_names()[i] != string("[Term]") ) continue;

	// Compose the node information in a temporary node
	bool obsolete = false;
	GONode tempNode;
	GOGraph* go_namespace = NULL;

	// The first tag-value pair must contain the id
	if( stanza_tags()[i][0] != string( "id" ) )
	  {
	    cout<<"Missing id. Skipping the stanza."<<endl;
	    continue;
	  }
	else
	  tempNode.id = stanza_vals()[i][0];

	// Traverse the tag-value pairs and compose node information
	for( unsigned int j = 0; j < stanza_tags()[i].size(); j++ )
	  {
	    // Ignore the obsolete terms
	    if( stanza_tags()[i][j] == string( "is_obsolete" ) &&
		stanza_vals()[i][j] == string( "true" ) )
	      { 
		obsolete = true;
		break;
	      }

	    // Parse the name
	    if( stanza_tags()[i][j] == string( "name" ) )
	      tempNode.name = stanza_vals()[i][j];

	    // Parse the namespace
	    if( stanza_tags()[i][j] == string( "namespace" ) )
	      {
		if( stanza_vals()[i][j] == string( "molecular_function" ) )
		  go_namespace = &go_mf;
		else if( stanza_vals()[i][j] == string( "biological_process" ) )
		  go_namespace = &go_bp;
		else if( stanza_vals()[i][j] == string( "cellular_component" ) )
		  go_namespace = &go_cc;
		else
		  cout<<"Unknown namespace for ID "<<tempNode.id<<endl;
	      }

	    // Parse the hierarchical relationships
	    if( stanza_tags()[i][j] == string( "is_a" ) )
	      {
		string str = stanza_vals()[i][j];

		// Remove the comment
		string::size_type k = str.find( "!" );
		if( k != string::npos ) str = str.substr( 0, k-1 );

		// Store the edge
		std::pair< string, string > edge( tempNode.id, str );
		edges.insert( edge );
	      }

	    if( pof == true && stanza_tags()[i][j] == string( "relationship" ) )
	      {
		string str = stanza_vals()[i][j];
	      
		// Remove the comment
		string::size_type k = str.find( "!" );
		if( k != string::npos ) str = str.substr( 0, k-1 );

		// Determine if this is a part_of relationship
		string::size_type pof_loc = str.find( "part_of" );
		if( pof_loc != string::npos )
		  {
		    // Store the edge
		    str = str.substr( pof_loc + 8 );
		    std::pair< string, string > edge( tempNode.id, str );
		    edges.insert( edge );
		  }
	      }
	  }

	if( obsolete == true ) continue;
	if( go_namespace == NULL ) continue;

	// Add the new vertex and fill in its property values
	GOvd vd = boost::add_vertex( *go_namespace );
	(*go_namespace)[vd] = tempNode;
      
	// Store the descriptor into the lookup table
	std::pair< GOvd, GOGraph* > sat_info( vd, go_namespace );
	std::pair< string, std::pair<GOvd, GOGraph*> > new_desc( tempNode.id, sat_info );
	vdesc.insert( new_desc );
      }

    // Add edges
    for( std::multimap< string, string >::iterator iter = edges.begin();
	 iter != edges.end(); iter++ )
      {
	// Find the vertex descriptors
	typedef map< string, std::pair<GOvd, GOGraph*> >::iterator GOvdi;
	GOvdi vdi1 = vdesc.find( iter->first );
	if( vdi1 == vdesc.end() )
	  {
	    cout<<"Edge construction: Failed to find "<<iter->first<<endl;
	    return false;
	  }

	GOvdi vdi2 = vdesc.find( iter->second );
	if( vdi2 == vdesc.end() )
	  {
	    cout<<"Edge construction: Failed to find "<<iter->second<<endl;
	    return false;
	  }

	// Make sure the nodes belong to the same ontology
	if( vdi1->second.second != vdi2->second.second )
	  {
	    cout<<"Edge construction: ";
	    cout<<iter->first<<" and "<<iter->second<<" are in different ontologies"<<endl;
	    return false;
	  }
      
	// Retrieve the vertex descriptors
	GOvd vd1 = vdi1->second.first;
	GOvd vd2 = vdi2->second.first;

	// Add the edge
	boost::add_edge( vd1, vd2, *(vdi1->second.second) );
      }

    // Clear up the temporary storage
    header_tags().clear();
    header_vals().clear();
    stanza_names().clear();
    stanza_tags().clear();
    stanza_vals().clear();

    return true;
  }

  // Outputs the container in GraphViz .dot format to an output stream
  void CGOContainer::toDot( std::ostream& os, ONTOLOGY_INDEX oi ) const
  {
    // Retrieve the ontology
    const GOGraph* g = index_ontology( oi );
    if( g == NULL ) return;

    // Write to the stream
    go_label_writer glw( *g );
    boost::write_graphviz( os, *g, glw );
  }

  // Outputs the container in GraphViz .dot format
  void CGOContainer::toDot( const char* filename, ONTOLOGY_INDEX oi ) const
  {
    std::ofstream ofs( filename );
    toDot( ofs, oi );
    ofs.close();
  }

  // Returns the number of nodes in the ontology
  unsigned int CGOContainer::size( ONTOLOGY_INDEX oi ) const
  {
    unsigned int res = 0;

    if( (oi & GO_MF) == GO_MF ) res += boost::num_vertices( go_mf );
    if( (oi & GO_BP) == GO_BP ) res += boost::num_vertices( go_bp );
    if( (oi & GO_CC) == GO_CC ) res += boost::num_vertices( go_cc );

    return res;
  }

  // Returns the ontology associated with the provided index
  const CGOContainer::GOGraph* CGOContainer::index_ontology( ONTOLOGY_INDEX oi ) const
  {
    if( (oi & GO_MF) == GO_MF ) return &go_mf;
    else if( (oi & GO_BP) == GO_BP ) return &go_bp;
    else if( (oi & GO_CC) == GO_CC ) return &go_cc;
    else return NULL;
  }

  // Reverse ontology indexing
  ONTOLOGY_INDEX CGOContainer::rev_index_ontology( const GOGraph* g ) const
  {
    if( g == &go_mf ) return GO_MF;
    else if( g == &go_bp ) return GO_BP;
    else if( g == &go_cc ) return GO_CC;
    else return GO_NONE;
  }

  // Returns the ontology that the ID belongs to
  ONTOLOGY_INDEX CGOContainer::find_ontology( string id ) const
  {
    GOcvdi vdi = vdesc.find( id );
    if( vdi == vdesc.end() ) return GO_NONE;
    else return rev_index_ontology( vdi->second.second );
  }

  // Takes a GO ID and returns the GO ID of its parents
  void CGOContainer::getParents( const string& id, vector< string >& res ) const
  {
    res.clear();

    // Retrieve the vertex descriptor and the corresponding graph
    GOcvdi iter = vdesc.find( id );
    if( iter == vdesc.end() ) return;
    GOvd vd = iter->second.first;
    GOGraph* g = iter->second.second;

    // Retrieve the outgoing edges
    std::pair< GOoei, GOoei > edges = boost::out_edges( vd, *g );

    // Traverse the outgoing edges
    for( GOoei e_iter = edges.first; e_iter != edges.second; e_iter++ )
      {
	// Retrieve the node on the target side of the edge
	GOvd t = boost::target( *e_iter, *g );

	// Store the parent ID
	res.push_back( (*g)[t].id );
      }
  }


  // Takes a GO ID and returns the full path to the root (including the node itself)
  void CGOContainer::getFullPath( const string& id, set< string >& res ) const
  {
    // Retrieve the vertex descriptor
    GOcvdi iter = vdesc.find( id );
    if( iter == vdesc.end() ) return;
    GOvd vd = iter->second.first;
    GOGraph* g = iter->second.second;

    // Add self
    res.insert( id );

    // Only consider vertices with outgoing edges
    if( boost::out_degree( vd, *(iter->second.second) ) < 1 ) return;

    // Retrieve the outgoing edges
    std::pair< GOoei, GOoei > edges = boost::out_edges( vd, *g );

    // Traverse the outgoing edges
    for( GOoei e_iter = edges.first; e_iter != edges.second; e_iter++ )
      {
	// Retrieve the node on the target side of the edge
	GOvd t = boost::target( *e_iter, *g );

	// Recurse
	getFullPath( (*g)[t].id, res );
      }
  }

  // Same thing but takes a vector of IDs, and returns the union of pathss
  void CGOContainer::getFullPaths( const vector< string >& ids, set< string >& res ) const
  {
    for( unsigned int i = 0; i < ids.size(); i++ )
      getFullPath( ids[i], res );
  }

  // Projects a set of paths onto the current ontology
  void CGOContainer::projectPaths( const set< string >& input, set< string >& output ) const
  {
    for( set< string >::const_iterator iter = input.begin();
	 iter != input.end(); iter++ )
      {
	if( this->find_ontology( *iter ) == GO_NONE ) continue;
	else output.insert( *iter );
      }
  }

  // Takes a vector of IDs, and returns the leaf nodes among
  void CGOContainer::getLeafs( const vector< string >& ids, vector< string >& leafs ) const
  {
    // Generate the parents set
    set< string > par;
    for( unsigned int i = 0; i < ids.size(); i++ )
      {
	vector< string > parents;
	getParents( ids[i], parents );
	for( unsigned int j = 0; j < parents.size(); j++ )
	  par.insert( parents[j] );
      }

    // Traverse the ids a second time and identify non-parent nodes; these are leafs
    for( unsigned int i = 0; i < ids.size(); i++ )
      {
	if( par.find( ids[i] ) == par.end() )
	  leafs.push_back( ids[i] );
      }
  }

  // Counts the number of nodes that two sets of paths share in common
  int CGOContainer::nPathNodesShared( const set< string >& s1,
				      const set< string >& s2 ) const
  {
    int res = 0;

    // Count the number of matching entries
    set< string >::const_iterator iter1 = s1.begin();
    set< string >::const_iterator iter2 = s2.begin();
  
    while( iter1 != s1.end() && iter2 != s2.end() )
      {
	// Matching labels
	if( *iter1 == *iter2 )
	  {
	    res++;
	    iter1++;
	    iter2++;
	  }

	// Iter 1 lagging behind
	else if( *iter1 < *iter2 ) iter1++;

	// Iter 2 lagging behind
	else iter2++;
      }

    return res;
  }
}
