
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SDC.H"
#include "gauss_radau.H"

namespace sdc
{
  SDC::SDC(
    std::shared_ptr<SDCSolver> solver,
    int nbNodes
    )
    :
    solver( solver ),
    nbNodes( nbNodes ),
    nodes(),
    smat(),
    qmat(),
    dsdc()
  {
    assert( solver );
    assert( nbNodes > 1 );
    assert( nbNodes < 14 );

    quadrature::rules( nbNodes, nodes, smat, qmat );

    dsdc.resize( nodes.rows() - 1 );

    for ( int i = 0; i < dsdc.rows(); i++ )
      dsdc( i ) = nodes( i + 1 ) - nodes( i );
  }

  SDC::~SDC()
  {}
}
