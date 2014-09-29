/// \file params.cpp
/// \brief Implementation of the parameter handling
/// \author Artem Sokolov

#include "params.h"

#include <fstream>
#include <iostream>

#include <boost/spirit/include/classic.hpp>
#include <boost/spirit/include/classic_file_iterator.hpp>
#include <boost/spirit/include/classic_increment_actor.hpp>
#include <boost/spirit/include/classic_clear_actor.hpp>

using std::cout;
using std::endl;

namespace spirit = boost::spirit::classic;

// Constructor that creates the object with default values
CStrutParams::CStrutParams()
  : myopts( "Supported options" ),
    exp_type_( "test" ),
    alg_choice_( "random" ),
    ioker_choice_( "prod" ),
    log_name_( "" )
{
  // Declare the available options
  myopts.add_options()
    ("exp_type", po::value< string >(&exp_type_), "Experiment type (one of test, ps, var)")
    ("alg_choice", po::value< string >(&alg_choice_), "Choice of the algorithm (one of prcp, 1svmm, 1svms, nsvmm, nsvms)")
    ("ioker_choice", po::value< string >(&ioker_choice_), "Choice of the IO joint kernel (one of prod, poly, polyh)")
    ("folds", po::value< string >(), "The folds withheld for testing during cross-validation" )
    ("alg_params", po::value< string >(), "Algorithm parameters" )
    ("log_name", po::value< string >(&log_name_), "(Part of) A filename for additional logging" );
}

void CStrutParams::load( const char* filename )
{
  // Read in and store the option variable map
  po::variables_map vm;
  std::ifstream ifs( filename );
  store( parse_config_file(ifs, myopts), vm );
  notify( vm );
  ifs.close();

  string str;

  // Parse the fold values
  if( vm.count("folds") > 0 )
    {
      str = vm["folds"].as<string>();
      if( spirit::parse( str.c_str(), *(spirit::uint_p[spirit::push_back_a(folds_)] >> *(spirit::blank_p)) ).full == false )
	throw boost::program_options::invalid_option_value("Unparsable fold values");
    }

  // Parse the algorithm parameters
  if( vm.count("alg_params") > 0 )
    {
      str = vm["alg_params"].as<string>();
      if( spirit::parse( str.c_str(), *(spirit::real_p[spirit::push_back_a(alg_params_)] >> *(spirit::blank_p)) ).full == false )
	throw boost::program_options::invalid_option_value("Unparsable fold values");
    }

  // Process any additional options
  processArgs( vm );
}

void CStrutParams::display() const
{
  cout<<"Experiment type      : "<<exp_type()<<endl;
  cout<<"Algorithm choice     : "<<alg_choice()<<endl;
  cout<<"IO kernel choice     : "<<ioker_choice()<<endl;
  cout<<"Folds                : "; 
  {
    vector< unsigned int > v = folds();
    std::copy( v.begin(), v.end(),
	       std::ostream_iterator<unsigned int>( cout, " " ) ); cout<<endl;
  }
  cout<<"Algorithm parameters : ";
  {
    vector< double > v = alg_params();
    std::copy( v.begin(), v.end(),
	       std::ostream_iterator<double>( cout, " " ) ); cout<<endl;
  }
  cout<<"Log name             : "<<log_name()<<endl;
}

////////////////////////////////////////////////////////////////////
CGOStrutParams::CGOStrutParams()
  : CStrutParams(),
    obo_location_("")
{
  ontology_ = GO::GO_MF;
  myopts.add_options()
    ("obo_location", po::value< string >(&obo_location_), "Full path to the .obo ontology")
    ("go_namespace", po::value< string >(), "Gene Ontology namespace associated with the experiment" );
}

void CGOStrutParams::processArgs( const boost::program_options::variables_map& vmap )
{
  ontology_ = GO::GO_NONE;

  if( vmap.count("go_namespace") == 0 ) return;
  string str = vmap["go_namespace"].as<string>();
  if( str.find("mf") != string::npos ) ontology_ |= GO::GO_MF;
  if( str.find("bp") != string::npos ) ontology_ |= GO::GO_BP;
  if( str.find("cc") != string::npos ) ontology_ |= GO::GO_CC;
}

void CGOStrutParams::display() const
{
  CStrutParams::display();
  cout<<"GO Ontologies        :";
  if( GO::hasMF( ontology() ) ) cout<<" MF";
  if( GO::hasBP( ontology() ) ) cout<<" BP";
  if( GO::hasCC( ontology() ) ) cout<<" CC";
  cout<<endl;
  cout<<"Path to .obo file    : "<<obo_location_<<endl;
}

