
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDC.H"

namespace sdc
{

  SDC::SDC( std::shared_ptr<SDCSolver> solver, int nbNodes )
  :
  solver( solver ),
  nbNodes( nbNodes )
  {
    assert( solver );
    assert( nbNodes > 1 );
    assert( nbNodes < 14 );
  }

  SDC::~SDC()
  {
    
  }

}
