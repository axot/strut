// -*-c++-*-
/// \file blastout.h
/// \brief headers for CBLASTOutput, storage container
/// \author Artem Sokolov

#ifndef BLASTOUT_H__INCLUDED
#define BLASTOUT_H__INCLUDED

#include "types.h"

namespace GO
{
  /// E-values below this threshold are considered significant
  const double e_val_threshold = 1e-6;

  /// Parses and stores tabular output from blastpgp
  class CBLASTOutput
  {
  public:
    /// Constructor
    CBLASTOutput( const string& filename );

  public:
    /// Stores entries of BLAST output
    class SBOEntry
    {
    public:
      /// Constructor
      SBOEntry( std::string raw_line );

    public:

      /// Query ID
      std::string query_id;

      /// Subject ID
      std::string subject_id;

      /// Percent of identity
      double percent_identity;

      /// Alignment length
      int alignment_length;

      /// Number of mismatches (not including gaps)
      int mismatches;

      /// Number of gap openings
      int gaps;

      /// Start of alignment in query
      int q_start;

      /// End of alignment in query
      int q_end;

      /// Start of alignment in subject
      int s_start;

      /// End of alignment in subject
      int s_end;

      /// Expected value
      double e_value;

      /// Bit score
      double bit_score;

    public:
      /// Displays the entry to an arbitrary output stream
      void display( std::ostream& os = std::cout ) const;

    public:
      /// Order by e-value
      bool operator< ( const SBOEntry& other ) const { return this->e_value < other.e_value; }
    };

  public:
    /// Iterator over the entries
    typedef std::map< std::string, std::vector< SBOEntry > >::const_iterator emap_const_iter_t;

  private:
    /// Maps a query ID to its respective set of entries
    std::map< std::string, std::vector< SBOEntry > > entries;

  public:
    /// Displays the entire BLAST output table to an arbitrary output stream
    void display( std::ostream& os = std::cout ) const;

    /// Returns the number of entries
    unsigned int size() const { return entries.size(); }

    /// Returns an iterator to the beginning of the table
    emap_const_iter_t begin() const { return entries.begin(); }

    /// Returns an iterator to the end of the table
    emap_const_iter_t end() const { return entries.end(); }

    /// Returns an iterator to the entry that matches the provided key (query ID)
    emap_const_iter_t find( const std::string& key ) const { return entries.find( key ); }

    /// Returns the hits associated with the requested key
    std::vector< SBOEntry > getHits( const std::string& key ) const;

    /// Returns the names of all entries
    std::vector< string > getNames() const;
  
    /// Returns true if two proteins are within e_thresh of each other
    bool proximityEVal( const string& p1, const string& p2, double e_thresh ) const;

    /// Returns true if two proteins are within pi_thresh of each other
    bool proximityPIden( const string& p1, const string& p2, double pi_thresh ) const;

  };
}

#endif
