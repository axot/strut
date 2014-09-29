// -*-c++-*-
/// \file params.h
/// \brief Parameter handling
/// \author Artem Sokolov

#ifndef PARAMS_H__INCLUDED
#define PARAMS_H__INCLUDED

#include <iostream>

#include "go-container.h"
#include "types.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

/// Generic parameters for the strut framework
class CStrutParams
{
public:
  /// Constructor that creates the object with default values
  CStrutParams();

  /// Destructor
  virtual ~CStrutParams() {}

  /// Parses parameters from a file
  void load( const char* filename );

private:
  /// Process any additional options by the child classes
  virtual void processArgs( const boost::program_options::variables_map& vmap ) {}

protected:
  po::options_description myopts;

private:
  /// The type of experiment
  string exp_type_;

  /// Choice of the algorithm
  string alg_choice_;

  /// Choice of the IO joint kernel
  string ioker_choice_;

  /// The fold(s) withheld for testing
  vector< unsigned int > folds_;

  /// The algorithm parameters
  vector< double > alg_params_;

  /// The filename for additional logging
  string log_name_;

public:

  const string exp_type() const
  {return exp_type_;}

  const string alg_choice() const
  {return alg_choice_;}

  const string ioker_choice() const
  {return ioker_choice_;}

  const vector< unsigned int > folds() const
  {return folds_;}

  const vector< double > alg_params() const
  {return alg_params_;}

  const string log_name() const
  {return log_name_;}

public:
  virtual void display() const;

  virtual void displayHelp() const
  {cout<<myopts<<endl;}
  
};

/// Specialized parameters for Gene Ontology experiments
class CGOStrutParams : public CStrutParams
{
public:
  /// Constructor
  explicit CGOStrutParams();

  /// Destructor
  virtual ~CGOStrutParams() {}
  
private:
  /// Process any additional options by the child classes
  virtual void processArgs( const boost::program_options::variables_map& vmap );

private:
  string obo_location_;
  GO::ONTOLOGY_INDEX ontology_;

public:
  const std::string obo_location() const
  {return obo_location_;}

  const GO::ONTOLOGY_INDEX ontology() const
  {return ontology_;}

public:
  virtual void display() const;
};

#endif

