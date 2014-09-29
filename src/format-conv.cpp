// format_conv.cpp - generates a sparse dataset from a BLAST or annotation file
// 
// by Artem Sokolov

#include "go-annotation.h"
#include "blastout.h"
#include "sample.h"

#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <boost/program_options.hpp>

using std::cout;
using std::endl;
using std::set;
using std::string;
using std::vector;
using std::shared_ptr;

using namespace GO;
namespace po = boost::program_options;

int main( int argc, char* argv[] )
{
  string format = "";
  string obo_location = "";
  string input_filename = "";
  string output_filename = "";

  // Compose the allowed program arguments
  po::options_description opts( "Supported options" );
  opts.add_options()
    ("input-format,f", po::value< string >( &format ), "Format of the input file. Must be goa or blast" )
    ("obo-location,b", po::value< string >( &obo_location ), "Ontology file to be used with goa files" )
    ("input,i", po::value< string >( &input_filename ), "Name of the input file" )
    ("output,o", po::value< string >( &output_filename ), "Name of the output file" );

  // Display usage as needed
  if( argc < 2 )
    {
      cout<<"Usage: "<<argv[0]<<" options"<<endl;
      cout<<opts;
      return -1;
    }

  // Parse the arguments
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, opts), vm);
  po::notify(vm);

  // Display the parsed arguments
  cout<<"Wokring with the following arguments: "<<endl;
  cout<<"Input  file: "<<input_filename<<endl;
  cout<<"Output file: "<<output_filename<<endl;
  cout<<".obo   file: "<<obo_location<<endl;
  cout<<"Format     : "<<format<<endl;
  
  // Constants
  const ONTOLOGY_INDEX myfilter = GO_MF | GO_BP | GO_CC;
  const double e_val_lower_thresh = 1e-10;
  const double e_val_upper_thresh = 50.0;

  // Dataset object to store the data
  CDataSet<CSparseSample> ds;

  if( format == "goa" )
    {
      // Pre-load the .obo file
      GO::CGOContainer goGraph( obo_location );

      // Parse a GO annotation file
      cout<<"Loading "<<input_filename<<"... "; cout.flush();
      CGOACollection goa( input_filename.c_str() );
      cout<<goa.size()<<" annotations parsed"<<endl;

      // Retrieve the names of proteins that have legit annotations
      vector< string > pids = goa.getAnnotatedProteins( myfilter );
      cout<<pids.size()<<" proteins have legit annotations"<<endl;

      // Make a sparse dataset out of the annotations
      cout<<"Generating the dataset... "; cout.flush();
      makeSparseDataset( goa, pids, ds, goGraph, myfilter );
      cout<<" generated "<<ds.size()<<" samples"<<endl;
    }
  else if( format == "blast" )
    {
      // Parse the BLAST file
      cout<<"Loading "<<input_filename<<"... "; cout.flush();
      CBLASTOutput blout( input_filename.c_str() );
      cout<<blout.size()<<" entries parsed"<<endl;
      
      // Make a sparse dataset out of the BLAST hits
      cout<<"Generating the dataset... "; cout.flush();
      makeSparseDataset( blout, ds, e_val_lower_thresh, e_val_upper_thresh );
      cout<<" generated "<<ds.size()<<" samples"<<endl;
    }
  else
    throw std::runtime_error( "Unknown source format. Options: goa, blast" );

  // Save it to a file in sparse format
  ds.save( output_filename );

  return 0;
}
