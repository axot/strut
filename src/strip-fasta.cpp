// strip-fasta.cpp - strips a fasta file to contain only those proteins that have 
//	matching annotations in the provided anootation file
//
// by Artem Sokolov

#include "types.h"
#include "go-annotation.h"

#include <iostream>

using std::logic_error;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

int main( int argc, char* argv[] )
{
  // Verify the arguments
  if( argc < 5 )
    {
      cout<<"Usage: "<<argv[0]<<" <input FASTA file> <annotations file> <output FASTA file> <CC|BP|MF>"<<endl;
      return -1;
    }

  // Parse the filter
  unsigned int iFil = boost::lexical_cast< unsigned int >( argv[4] );
  if( iFil > 7 )
    throw logic_error( "GO filter must be provided as a binary mask for CC|BP|MF and must be an integer between 0 and 7" );
  GO::ONTOLOGY_INDEX filter( iFil );

  // Display the arguments
  cout<<"Input  file: "<<argv[1]<<endl;
  cout<<"Output file: "<<argv[3]<<endl;
  cout<<"Annots file: "<<argv[2]<<endl;
  cout<<"GO   filter:";
  if( GO::hasMF( filter ) ) cout<<" MF";
  if( GO::hasBP( filter ) ) cout<<" BP";
  if( GO::hasCC( filter ) ) cout<<" CC";
  cout<<endl;

  // Read in the annotations file
  cout<<"Loading the GO annotations... "; cout.flush();
  GO::CGOACollection goa( argv[2] );
  cout<<"parsed "<<goa.size()<<" annotations"<<endl;

  // Open the input and output files
  ifstream ifs( argv[1] );
  ofstream ofs( argv[3] );

  // Traverse the input file
  bool copyMode = false;
  string str;
  int iTotal = 0;
  int iFilt = 0;
  for( std::getline( ifs, str ); ifs.fail() == false; std::getline( ifs, str ) )
    {
      // Is this a header?
      if( str[0] == '>' )
	{
	  if( (++iTotal % 1000) == 0 ) {cout<<"."; cout.flush();}
	  copyMode = false;

	  // Determine if the header has any significant annotations
	  string name = str.substr( 1 );
	  string id = goa.getObjectID( name, false );
	  if( goa.hasGOIDs( id, filter ) ) {copyMode = true; ++iFilt;}
	}

      if( copyMode == true )
	ofs<<str<<endl;
    }
  cout<<endl;
  cout<<"Filtered "<<iTotal<<" sequences down to "<<iFilt<<" sequences"<<endl;

  // Close the files
  ifs.close();
  ofs.close();

  return 0;
}
