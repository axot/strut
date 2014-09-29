// cafa-post.cpp - post-processing of CAFA predictions
//
// by Artem Sokolov

#include "types.h"
#include "go-annotation.h"
#include "sample.h"
#include "parsers.h"

#include <fstream>

#include <boost/algorithm/string.hpp>

typedef struct
{
  string name;
  double pIdent;
  unsigned int mLength;
} id_entry;

typedef std::unordered_map< string, id_entry > idmap_t;

// Loads the ID map
idmap_t loadIDMap( const string& filename )
{
  idmap_t res;

  // Open the file and traverse it line-by-line
  std::ifstream ifs( filename.c_str() );
  string s;
  for( std::getline( ifs, s ); ifs.fail() == false; std::getline( ifs, s ) )
    {
      // Tokenize the string
      vector< string > toks;
      boost::split( toks, s, boost::is_space() );
      if( toks.size() < 4 ) throw std::logic_error( "Invalid format of the id match file" );

      // Parse percent-identy and match-length
      double pi = boost::lexical_cast<double>( toks[2] );
      double ml = boost::lexical_cast<unsigned int>( toks[3] );

      // Handle existing entry
      if( res.find( toks[0] ) != res.end() )
	{
	  // Percent identity needs to be higher for replacement
	  if( pi < res[toks[0]].pIdent ) continue;

	  // Match length needs to be higher for replacement
	  if( ml < res[toks[0]].mLength ) continue;
	}

      // Add / replace the entry
      id_entry e = { toks[1], pi, ml };
      res[toks[0]] = e;
    }

  return res;
}

// Maps a set of IDs in a dataset using the provided ID map
/*set< string > mapIDs( const CDataSet<CSparseSample>& ds, const string& fn )
{

  // Traverse the dataset IDs and use the mapping
  set< string > res;
  for( unsigned int i = 0; i < ds.size(); ++i )
    {
      string s = ds.i2s( i );
      idmap_t::const_iterator iter = m.find( s );
      if( iter == m.end() ) continue;
      res.insert( iter->second.name );
    }
  return res;
  }*/

int main( int argc, char* argv[] )
{
  // Display usage
  if( argc < 7 )
    {
      cout<<"Usage: "<<argv[0]<<" <prediction file> <id-match file> <annotation file> <.obo GO ontology file> <output file> <ontology> [ppi-based prediction file]"<<endl;
      return -1;
    }

  shared_ptr< CFeatMap > pfm( new CFeatMap );
  CDataSet<CSparseSample> dsPred;

  // Load the prediction data
  if( argc == 8 )
    {
      cout<<"Loading the ppi-based predictions... "; cout.flush();
      parseSparseFile( argv[7], dsPred, pfm, ',', ',', '=' );
      cout<<"parsed "<<dsPred.size()<<" samples, "<<nFeats( dsPred )<<" features"<<endl;
    }

  // Load the non-ppi-based prediction data (overwrite == false)
  cout<<"Loading the predictions... "; cout.flush();
  parseSparseFile( argv[1], dsPred, pfm, ',', ',', '=' );
  cout<<"data now has "<<dsPred.size()<<" samples, "<<nFeats( dsPred )<<" features"<<endl;

  // Find a set of ids matched to the test data
  // Load the id map
  cout<<"Loading the id map... "; cout.flush();
  idmap_t m = loadIDMap( argv[2] );
  cout<<"loaded "<<m.size()<<" entries"<<endl;

  // Loading the Ontology file
  GO::CGOContainer goGraph( argv[4] );

  // Identify the ontology
  string sOntology( argv[6] );
  GO::ONTOLOGY_INDEX goFilter;
  if( sOntology.find( "mf" ) != sOntology.npos )
    {
      cout<<"Focusing on molecular function"<<endl;
      goFilter = GO::GO_MF;
    }
  else if( sOntology.find( "bp" ) != sOntology.npos )
    {
      cout<<"Focusing on biological process"<<endl;
      goFilter = GO::GO_BP;
    }
  else throw std::logic_error( "Ontology not supported" );

  // Load the GO annotation file
  cout<<"Loading the GO annotations... "; cout.flush();
  GO::CGOACollection goa( argv[3] );
  cout<<"loaded "<<goa.size()<<" annotations"<<endl;

  std::ofstream ofs( argv[5] );

  ofs<<"AUTHOR GOstruct"<<endl;
  if( argc == 8 )
    {
      ofs<<"MODEL 1"<<endl;
      ofs<<"KEYWORDS sequence alignments, sequence properties, protein interactions, machine learning based method."<<endl;
    }
  else
    {
      ofs<<"MODEL 2"<<endl;
      ofs<<"KEYWORDS sequence alignments, sequence properties, machine learning based method."<<endl;
    }

  // Traverse the predictions
  for( unsigned int i = 0; i < dsPred.size(); ++i )
    {
      set< string > annots;
      
      // Retrieve the annotations associated with this prediction
      shared_ptr< CSparseSample const > p = dsPred.getSample( i );
      for( unsigned int j = 0; j < pfm->nFeats(); ++j )
	{
	  string f = pfm->i2f( j );
	  if( p->getValue( f ) != 0.0 )
	    annots.insert( f );
	}
  
      if( annots.empty() ) continue;

      // Find a match for the protein
      string s = dsPred.i2s( i );
      idmap_t::const_iterator siter = m.find( s );
      if( siter != m.end() )
	{
	  // Retrieve the set of GO IDs for the protein and project them onto the slim ontology
	  set< string > annot_full;
	  vector< string > a = goa.getGOIDs( siter->second.name, goFilter );

	  goGraph.getFullPaths( a, annot_full );
	  for( set< string >::iterator iter = annot_full.begin();
	       iter != annot_full.end(); ++iter )
	    annots.insert( *iter );
	}


      // Output the predictions
      for( set<string>::iterator iter = annots.begin();
	   iter != annots.end(); ++iter )
	ofs<<s<<" "<<*iter<<" 1.00"<<endl;
    }

  ofs<<"END"<<endl;
  
  return 0;
}
