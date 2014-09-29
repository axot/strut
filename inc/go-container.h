// -*-c++-*-
/// \file go-container.h
/// \brief headers for CGOContainer
/// \author Artem Sokolov

#ifndef GO_CONTAINER_H__INCLUDED
#define GO_CONTAINER_H__INCLUDED

#define BOOST_NO_HASH
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>

#include "types.h"

namespace GO
{
  typedef int ONTOLOGY_INDEX;

  extern const ONTOLOGY_INDEX GO_NONE;
  extern const ONTOLOGY_INDEX GO_MF;	// Molecular function
  extern const ONTOLOGY_INDEX GO_BP;	// Biological process
  extern const ONTOLOGY_INDEX GO_CC;	// Cellular component

  bool hasMF( ONTOLOGY_INDEX i );
  bool hasBP( ONTOLOGY_INDEX i );
  bool hasCC( ONTOLOGY_INDEX i );

  /// Gene ontology container
  class CGOContainer
  {
  private:

    /// Graph node information
    struct GONode
    {
      /// The unique GO ID of the term
      string id;

      /// The term name
      string name;
    };

    /// The graph definition. The edges point from children to parents
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, GONode> GOGraph;

    /// The associated vertex descriptor
    typedef boost::graph_traits<GOGraph>::vertex_descriptor GOvd;

    /// The associated vertex iterator
    typedef boost::graph_traits<GOGraph>::vertex_iterator GOvi;

    /// The associated iterator over outgoing edges
    typedef boost::graph_traits<GOGraph>::out_edge_iterator GOoei;

  public:
    /// Constructor
    explicit CGOContainer( const string& filename ) {parseOBO( filename.c_str() );}

    /// Destructor
    ~CGOContainer() {}

  private:
    /// Molecular function ontology
    GOGraph go_mf;

    /// Biological process ontology
    GOGraph go_bp;
  
    /// Cellular component ontology
    GOGraph go_cc;

    /// Iterator over CGOContainer::vdesc
    typedef map< string, std::pair<GOvd, GOGraph*> >::iterator GOvdi;

    /// Const iterator over CGOContainer::vdesc
    typedef map< string, std::pair<GOvd, GOGraph*> >::const_iterator GOcvdi;

    /// Stores descriptors for efficient vertex search
    map< string, std::pair<GOvd, GOGraph*> > vdesc;

    /// A Label writer for .dot format
    class go_label_writer 
    {
    public:

      /// Constructor
      go_label_writer(const GOGraph& _g) : g(&_g) {}

      /// Writes out the label for a vertex or an edge
      template <class VertexOrEdge>
      void operator()(std::ostream& out, const VertexOrEdge& v) const {
	out << "[label=\"" << (*g)[v].id << "\"]";
      }

    private:
      /// The associated graph
      const GOGraph* g;
    };

  public:
    /// Parses an .obo file
    /** Set pof to true if "part_of" relationships should be counted as "is_a" for edge construction
	\return false in case of a failure, true otherwise
    */
    bool parseOBO( const char* filename, bool pof = true );

    /// Outputs the container in GraphViz .dot format to an arbitrary stream
    void toDot( std::ostream& os, ONTOLOGY_INDEX oi ) const;

    /// Outputs the container in GraphViz .dot format to a file
    void toDot( const char* filename, ONTOLOGY_INDEX oi ) const;

    /// Returns the number of nodes in the ontology
    unsigned int size( ONTOLOGY_INDEX oi ) const;

    /// Returns the ontology that the ID belongs to
    ONTOLOGY_INDEX find_ontology( string id ) const;

    /// Takes a GO ID and returns the GO ID of its parents
    void getParents( const string& id,
		     vector< string >& res ) const;

    /// Takes a GO ID and returns the full path to the root (including the node itself)
    void getFullPath( const string& id, set< string >& res ) const;

    /// Takes a vector of IDs, and returns the union of pathss from each ID to the root
    void getFullPaths( const vector< string >& ids, set< string >& res ) const;

    /// Takes a vector of IDs, and returns the leaf nodes among
    void getLeafs( const vector< string >& ids, vector< string >& leafs ) const;

    /// Projects a set of paths onto the current ontology
    void projectPaths( const set< string >& input, set< string >& output ) const;

    /// Counts the number of nodes that two sets of paths share in common
    int nPathNodesShared( const set< string >& s1,
			  const set< string >& s2 ) const;

  private:
    /// Returns the ontology associated with the provided index
    const GOGraph* index_ontology( ONTOLOGY_INDEX oi ) const;

    /// Returns the index associated with the ontology
    ONTOLOGY_INDEX rev_index_ontology( const GOGraph* g ) const;
  };
}

#endif

