// gostruct-dataprep.cpp - prepares the gostruct data
//
// by Artem Sokolov

#include <iostream>
#include <iterator>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

#include "opt.h"

int main( int argc, char* argv[] )
{
  vector< double > Q(4);
  Q[0] = 0.211146;

  vector< double > b(2);
  b[0] = -0.117217;

  vector< double > x = strongQuadraticOpt( Q, b, 100.0 );

  return 0;
}

