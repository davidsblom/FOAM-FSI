
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "SpaceMapping.H"

using namespace fsi;

SpaceMapping::SpaceMapping(
  shared_ptr<SurrogateModel> fineModel,
  shared_ptr<SurrogateModel> coarseModel,
  int maxIter,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex
  )
  :
  fineModel( fineModel ),
  coarseModel( coarseModel ),
  maxIter( maxIter ),
  nbReuse( nbReuse ),
  reuseInformationStartingFromTimeIndex( reuseInformationStartingFromTimeIndex ),
  timeIndex( 0 ),
  coarseResiduals(),
  fineResiduals(),
  coarseResidualsList(),
  fineResidualsList()
{
  assert( fineModel );
  assert( coarseModel );
  assert( maxIter > 0 );
  assert( nbReuse >= 0 );
  assert( reuseInformationStartingFromTimeIndex >= 0 );
};

void SpaceMapping::finalizeTimeStep()
{
  timeIndex++;
}

bool SpaceMapping::isConvergence(
  const fsi::vector & xk,
  const fsi::vector & xkprev,
  bool residualCriterium
  )
{
  if ( residualCriterium )
    return fineModel->isConvergence();

  return fineModel->isConvergence( xk, xkprev );
}

void SpaceMapping::iterationsConverged()
{
  // Save input/output information for next time step
  if ( nbReuse > 0 && fineResiduals.cols() >= 2 && timeIndex >= reuseInformationStartingFromTimeIndex )
  {
    assert( coarseResiduals.cols() == fineResiduals.cols() );

    coarseResidualsList.push_front( coarseResiduals );
    fineResidualsList.push_front( fineResiduals );
  }

  // Remove the last item from the residual list and solutions
  // list in order to ensure that at maximum nbReuse time steps
  // are included.
  while ( static_cast<int>( fineResidualsList.size() ) > nbReuse )
  {
    coarseResidualsList.pop_back();
    fineResidualsList.pop_back();
  }

  coarseResiduals.resize( coarseResiduals.rows(), 0 );
  fineResiduals.resize( fineResiduals.rows(), 0 );

  assert( coarseResidualsList.size() == fineResidualsList.size() );
  assert( static_cast<int>( coarseResidualsList.size() ) <= nbReuse );
}

void SpaceMapping::performPostProcessing(
  const fsi::vector & x0,
  fsi::vector & xk
  )
{
  fsi::vector y( x0.rows() );
  y.setZero();
  performPostProcessing( y, x0, xk, true );
}

void SpaceMapping::performPostProcessing(
  const fsi::vector & y,
  const fsi::vector & x0,
  fsi::vector & xk
  )
{
  performPostProcessing( y, x0, xk, false );
}
