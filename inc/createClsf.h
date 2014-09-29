// -*-c++-*-

#ifndef CREATECLSF_H__INCLUDED
#define CREATECLSF_H__INCLUDED

#include "random-clsf.h"
#include "perceptron.h"
#include "ssvm.h"
#include "nssvm.h"

/// Creates a classifier based on the parameters object
template< typename _I, typename _O >
shared_ptr< CClassifier<_I,_O> > createClassifier( const CStrutParams& params )
{
  shared_ptr< CClassifier<_I,_O> > pclsf;
  if( params.alg_choice() == string( "random" ) )
    pclsf.reset( new CRandomClassifier<_I,_O>() );
  else if( params.alg_choice() == string( "prcp" ) ||
	   params.alg_choice() == string( "prcp1") )
    {
      // Specify perceptron parameters
      SPerceptronParams pp;
      pp.bLossUpdate = (params.alg_choice() == string("prcp"));
      pp.margin = params.alg_params()[0];

      // Create the classifier
      pclsf.reset( new CPerceptron<_I,_O>( pp ) );
    }
  else
    {
      // Specify structural SVM parameters
      SSSVMParams svmp;
      svmp.Cn = params.alg_params()[0];
      svmp.eps = 0.01;
      svmp.nMaxQPSteps = 1000;
      svmp.fnPrefix = params.log_name();
      
      // Create the classifier
      if( params.alg_choice() == string( "1svmm" ) )
	pclsf.reset( new C1sSSVM_margin<_I, _O>( svmp ) );
      else if( params.alg_choice() == string( "1svms" ) )
	pclsf.reset( new C1sSSVM_slack<_I, _O>( svmp ) );
      else if( params.alg_choice() == string( "nsvmm" ) )
	pclsf.reset( new CnsSSVM<_I,_O,'m'>( svmp ) );
      else if( params.alg_choice() == string( "nsvms" ) )
	pclsf.reset( new CnsSSVM<_I,_O,'s'>( svmp ) );
      else
	throw std::logic_error( "Unrecognized classifier" );
    }

  return pclsf;
}

#endif

