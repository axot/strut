// -*-c++-*-
/// \file go-annotation.h
/// \brief headers for CGOAnnotation and CGOACollection
/// \author by Artem Sokolov

#ifndef GO_ANNOTATION_H__INCLUDED
#define GO_ANNOTATION_H__INCLUDED

#include "go-container.h"

namespace GO
{
  /// Returns true if the evidence code is safe
  bool isSafeEvidenceCode( const string& str );

  /// Gene Ontology Annotation
  class CGOAnnotation
  {
  public:
    /// Constructor
    explicit CGOAnnotation( const string& str );

  private:

    /// Field 1 (required): Database contributing the gene association
    string r_DB;

    /// Field 2 (required): A unique identifier for the item being annotated
    string r_ObjectID;

    /// Field 3 (required): A symbol to which ObjectID is matched
    string r_ObjectSym;

    /// Field 4 (optional): Flags that modify the interpretation
    string r_Qualifier;
 
    /// Field 5 (required): The GO annotation for ObjectID
    string r_GOID;
  
    /// Field 6 (required): Identification of the source
    string r_DBRef;
  
    /// Field 7 (required): Evidence code
    string r_ECode;
  
    /// Field 8 (optional): Additional identifiers
    string r_WorF;
  
    /// Field 9 (required): Which of the three ontology trees
    string r_Aspect;
  
    /// Field 10 (optional): Name of gene or gene product
    string r_ObjectName;
  
    /// Field 11 (optional): Synonymous identifiers
    string r_Synonym;
  
    /// Field 12 (required): Type of what is being annotated
    string r_ObjectType;
  
    /// Field 13 (required): Taxonomic identifiers
    string r_Taxon;
  
    /// Field 14 (required): Date when the annotation was made
    string r_Date;
  
    /// Field 15 (required): Database which made the annotation
    string r_AssignedBy;

    /// Field 16 (optional): Annotation extension
    string r_AnnotExt;

    /// Field 17 (optional): Gene Product Form ID
    string r_GenePFID;

  private:

  public:
  
    /// Displays the record to the provided output stream
    void display( std::ostream& os = std::cout ) const;

    /// Returns CGOAnnotation::r_ObjectID
    string getObjectID() const {return r_ObjectID;}

    /// Returns CGOAnnotation::r_ObjectSym
    string getObjectSymbol() const {return r_ObjectSym;}

    /// Returns CGOAnnotation::r_Qualifier
    string getQualifier() const {return r_Qualifier;}

    /// Returns CGOAnnotation::r_GOID
    string getGOID() const {return r_GOID;}

    /// Returns CGOAnnotation::r_ECode
    string getEvidenceCode() const {return r_ECode;}

    /// Returns CGOAnnotation::r_Synonym;
    string getSynonyms() const {return r_Synonym;}

    /// Retruns some combination of GO_MF, GO_BP and GO_CC
    ONTOLOGY_INDEX getAspect() const;

  };

  /// A collection of GO annotations for an organism
  class CGOACollection
  {
  private:
    /// Contains (ObjectID, (Object names, indices into annots)) entries
    /** i.e.,    (protein name, (list of protein synonyms, list of indices into annots)) entries
     */
    typedef map< string, std::pair< vector< string >, 
				    vector< unsigned int > > > goa_map_t;

    /// Const iterator for CGOACollection::goa_map_t
    typedef map< string, std::pair< vector< string >, 
				    vector< unsigned int > > >::const_iterator goa_map_citer_t;

  public:
    /// Constructor
    explicit CGOACollection( const string& filename );

  private:
    /// Raw annotations
    vector< CGOAnnotation > annots;

    /// Stores (protein, annotation info) entries
    goa_map_t goa_map;

    /// Names of all the proteins (including synonyms) in the collection
    set< string > prot_names;

  private:
    /// Post-processing of raw GO annotation values
    void postprocessing();

    /// Takes a name of a protein, searches through all synonyms
    /** \return An interator into CGOACollection::goa_map. Returns CGOACollection::goa_map.end() if no match found
     */
    goa_map_citer_t findProtein( const string& prot_name, bool bCaseSensitive = true ) const;

    /// Returns true if the j^th annotation is not a NOT, has a safe evidence code and satisfies the filter
    bool isAnnotationGood( unsigned int j, ONTOLOGY_INDEX filter ) const;

  public:

    /// Returns the i^th annotation
    CGOAnnotation getAnnot( unsigned int i ) const { return annots[i]; }

    /// Returns the number of annotations contained in the collection
    unsigned int size() const { return annots.size(); }

    /// Displays the i^{th} annotation
    void display( int i, std::ostream& os = std::cout ) const { annots[i].display( os ); }
  
    /// Displays the entire collection
    void display( std::ostream& os = std::cout ) const;

    /// Returns the object ID given its synonym or "" if there is no such synonym
    const string getObjectID( const string& syn, bool bCaseSensitive = true ) const;

    /// Returns a vector of GO IDs for a protein
    vector< string > getGOIDs( const string& prot_name,
			       ONTOLOGY_INDEX filter = (GO_MF | GO_BP | GO_CC),
			       bool bCaseSensitive = true ) const;

    /// Returns whether or not the protein contains GO IDs that match the filter
    bool hasGOIDs( const string& prot_name,
		   ONTOLOGY_INDEX filter = (GO_MF | GO_BP | GO_CC),
		   bool bCaseSensitive = true ) const;

    /// Returns all proteins that have a non-NOT annotation that falls into one of the provided categories
    vector< string > getAnnotatedProteins( ONTOLOGY_INDEX filter = (GO_MF | GO_BP | GO_CC) ) const;

    /// Appends the names of annotated proteins to res
    void getAnnotatedProteins( ONTOLOGY_INDEX filter, set< string >& res ) const;

  };
}

#endif
