// gostruct-dataprep.cpp - prepares the gostruct data
//
// by Artem Sokolov

#include "sample.h"
#include "blastout.h"
#include "go-annotation.h"
#include "io-dataset.h"

string annot_list_filename( "gostruct/foursp.annot_list" );
string blast_hits_filename( "/s/chopin/c/proj/protfun/data/BLAST/foursp/foursp.blast.gz" );

const GO::ONTOLOGY_INDEX goa_filter = GO::GO_MF;

// Loads a list of GOA collections from a file that specifies the collections' filenames
void loadGOACollections( string filename,
			 vector< std::shared_ptr< GO::CGOACollection > >& res )
{
  // Open the list of files
  std::ifstream ifs( filename.c_str() );
  string str;
  
  // Load in the datasets one at a time
  for( std::getline( ifs, str ); ifs.fail() == false; std::getline( ifs, str ) )
    {
      cout<<"Loading "<<str<<"..."; cout.flush();
      std::shared_ptr< GO::CGOACollection > pds( new GO::CGOACollection( str.c_str() ) );
      res.push_back( pds );
      cout<<pds->size()<<" annotations loaded"<<endl;
    }
}

// Returns the fold of each protein in a BLAST output file
// The folds are determined by the annotation files
// The value of -1 is used to denoted non-annotated proteins
void getBLASTfolds( shared_ptr< GO::CBLASTOutput const > blast_hits,
		    const std::vector< shared_ptr< GO::CGOACollection > >& annot_list,
		    simap_t& result )
{
  std::cout<<"Mapping BLAST hits to annotations..."<<std::endl;

  // Traverse the BLAST queries
  for( GO::CBLASTOutput::emap_const_iter_t iter = blast_hits->begin();
       iter != blast_hits->end(); iter++ )
    {
      // Retrieve the protein ID
      std::string protein_id = iter->first;

      // Determine which fold the protein in is by looking up its true annotation
      int iFold = -1;
      for( unsigned int i = 0; i < annot_list.size(); i++ )
	{
	  bool b_true_label = annot_list[i]->hasGOIDs( protein_id, goa_filter );
	  if( b_true_label == false ) continue;
	  else
	    {
	      iFold = i;
	      break;
	    }
	}

      result[ protein_id ] = iFold;
    }
}

// Returns the names of proteins that have significant BLAST hits and contain transferable annotations
// Outer index is into the fold (annot_list)
// Inner index is over the protein names in that GOACollection
void composeBLASTHitList( shared_ptr< GO::CBLASTOutput const > blast_hits,
			  const std::vector< shared_ptr< GO::CGOACollection > >& annot_list,
			  std::vector< std::set< std::string > >& res) 
{
  // Retrieve the number of folds
  unsigned int n = annot_list.size();
  res.resize( n );

  simap_t fold_map;
  getBLASTfolds( blast_hits, annot_list, fold_map );

  // Traverse the BLAST queries
  std::cout<<"Searching for transferable annotations in significant BLAST hits... "; std::cout.flush();
  for( GO::CBLASTOutput::emap_const_iter_t iter = blast_hits->begin();
       iter != blast_hits->end(); iter++ )
    {
      // Retrieve the protein ID and its fold
      std::string protein_id = iter->first;
      int iFold = fold_map[protein_id];

      // Skip non-annotated proteins
      if( iFold < 0 ) continue;

      // Retrieve the set of BLAST hits and sort them by e-value
      std::vector< GO::CBLASTOutput::SBOEntry > blast_hits = iter->second;
      std::sort( blast_hits.begin(), blast_hits.end() );

      // Traverse the BLAST hits until able to transfer the annotation
      bool b_pred_label = false;
      for( unsigned int j = 0; j < blast_hits.size(); j++ )
	{
	  // Make sure we haven't crossed the threshold yet for e-value yet
	  // No need to check the rest because the vector is sorted
	  if( blast_hits[j].e_value > GO::e_val_threshold ) break;

	  // Retrieve the match ID
	  std::string match_id = blast_hits[j].subject_id;

	  // Find the match's fold
	  int jFold = fold_map[match_id];

	  if( jFold >= 0 && jFold != iFold )
	    {
	      b_pred_label = true;

	      // Store the match to enforce symmetry
	      res[jFold].insert( match_id );
	      break;
	    }
	}
	      
      // Store the protein's name if it contains a transferable annotation
      if( b_pred_label == true )
	res[iFold].insert( protein_id );
    }

  // Display some statistics
  for( unsigned int i = 0; i < res.size(); i++ )
    std::cout<<res[i].size()<<" ";
  std::cout<<std::endl;
}

int main( int argc, char* argv[] )
{
  string str;

  // Verify the arguments
  if( argc < 2 )
    {
      cout<<"Usage: "<<argv[0]<<" <location of the .obo file>"<<endl;
      return -1;
    }

  // Pre-fetch ontology
  string obo_location( argv[1] );
  GO::CGOContainer goGraph( obo_location );

  // Load the annotations
  str = annot_list_filename;
  vector< shared_ptr< GO::CGOACollection > > annot_list;
  loadGOACollections( str, annot_list );

  // Load the blast hits
  str = blast_hits_filename;
  std::cout<<"Loading "<<str<<"..."; std::cout.flush();
  shared_ptr< GO::CBLASTOutput > blast_hits( new GO::CBLASTOutput( str.c_str() ) );
  std::cout<<blast_hits->size()<<" entries parsed"<<std::endl;

  // Retrieve proteins that have significant BLAST hits and transferable annotations
  vector< std::set< string > > good_prots;
  composeBLASTHitList( blast_hits, annot_list, good_prots );

  // Combine the names into a single set
  std::set< string > good_prots_all;
  for( unsigned int i = 0; i < good_prots.size(); i++ )
    good_prots_all.insert( good_prots[i].begin(), good_prots[i].end() );
  std::vector< string > good_prots_v( good_prots_all.begin(), good_prots_all.end() );

  // Generate the input space dataset
  std::cout<<"Generating a sparse dataset from BLAST hits... ";
  shared_ptr< CDataSet<CSparseSample> > pds_blast( new CDataSet<CSparseSample> );
  makeSparseDataset( *blast_hits, *pds_blast, 1e-10, 50.0 );
  pds_blast->subsample( good_prots_v );
  std::cout<<pds_blast->size()<<" samples"<<std::endl;
  pds_blast->save( "gostruct/input.sdat" );

  // Generate the output space datasets
  vector< shared_ptr< CDataSet<CSparseSample> > > pds_annot( good_prots.size() );
  vector< shared_ptr< CFeatMap > > pfms( good_prots.size() );
  for( unsigned int i = 0; i < good_prots.size(); i++ )
    {
      std::cout<<"Generating a sparse dataset from annotations... ";

      vector< string > myprots_v( good_prots[i].begin(), good_prots[i].end() );
	
      // Compose the output space dataset
      pds_annot[i].reset( new CDataSet<CSparseSample>() );
      pfms[i].reset( new CFeatMap );
      makeSparseDataset( *(annot_list[i]), myprots_v, *(pds_annot[i]), goGraph, goa_filter, pfms[i] );
      std::cout<<pds_annot[i]->size()<<" samples; "<<nFeats(pds_annot[i])<<" features"<<std::endl;
    }

  // Find all samples that have representation in both input and output spaces
  vector< string > sids1 = pds_blast->getSampleIDs();
  std::sort( sids1.begin(), sids1.end() );
  set< string > sids2;	// Sorted by construction
  for( unsigned int k = 0; k < good_prots.size(); k++ )
    for( unsigned int i = 0; i < pds_annot[k]->size(); i++ )
      sids2.insert( pds_annot[k]->i2s( i ) );
  set< string > sids;
  std::set_intersection( sids1.begin(), sids1.end(),
			 sids2.begin(), sids2.end(),
			 std::inserter( sids, sids.begin() ) );

  cout<<sids.size()<<" samples are represented in both input and output spaces"<<endl;

  // Count feature representation
  simap_t featRep;
  for( unsigned int k = 0; k < good_prots.size(); k++ )
    {
      // Traverse the output-space samples
      for( unsigned int i = 0; i < pds_annot[k]->size(); i++ )
	{
	  // Skip un-represented samples
	  string s = pds_annot[k]->i2s( i );
	  if( sids.find( s ) == sids.end() ) continue;

	  // Traverse the features
	  for( unsigned int j = 0; j < pfms[k]->nFeats(); j++ )
	    {
	      string f = pfms[k]->i2f( j );
	      if( pds_annot[k]->getSample(i)->getValue(j) != 0.0 )
		{
		  if( featRep.find(f) == featRep.end() )
		    featRep[f] = 1;
		  else
		    featRep[f] += 1;
		}
	    }
	}
    }

  // Find well-represented features
  int nThresh = 10;
  vector< string > goodFeats;
  for( simap_t::const_iterator iter = featRep.begin(); iter != featRep.end(); iter++ )
    if( iter->second >= nThresh ) goodFeats.push_back( iter->first );
  shared_ptr< CFeatMap > pfmGood( new CFeatMap( goodFeats ) );

  cout<<goodFeats.size()<<" features are well-represented"<<endl;

  // Subfeature the output space datasets
  for( unsigned int i = 0; i < good_prots.size(); i++ )
    {
      // Remove samples that consist of root-only annotations
      remap( *(pds_annot[i]), pfmGood );
      cropSamples( 2, *(pds_annot[i]) );

      // Save the dataset
      string filename = "gostruct/output" + boost::lexical_cast<string>(i) + ".sdat";
      pds_annot[i]->save( filename );
    }

  return 0;
}

