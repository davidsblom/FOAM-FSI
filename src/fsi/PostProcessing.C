
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PostProcessing.H"

using namespace fsi;

PostProcessing::PostProcessing(
  shared_ptr<MultiLevelFsiSolver> fsi,
  double initialRelaxation,
  int maxIter,
  int maxUsedIterations,
  int nbReuse,
  int reuseInformationStartingFromTimeIndex
  )
  :
  fsi( fsi ),
  initialRelaxation( initialRelaxation ),
  maxIter( maxIter ),
  maxUsedIterations( maxUsedIterations ),
  nbReuse( nbReuse ),
  timeIndex( 0 ),
  reuseInformationStartingFromTimeIndex( reuseInformationStartingFromTimeIndex ),
  residuals(),
  sols(),
  residualsList(),
  solsList(),
  residualsTimeList(),
  solsTimeList()
{
  assert( fsi );
  assert( maxIter > 0 );
  assert( initialRelaxation <= 1 );
  assert( initialRelaxation > 0 );
  assert( nbReuse >= 0 );
  assert( fsi->N > 0 );
  assert( !fsi->fluid->init );
  assert( !fsi->solid->init );
  assert( maxUsedIterations > 0 );
}

void PostProcessing::finalizeTimeStep()
{
  // Save input/output information for next time step
  if ( nbReuse > 0 && residualsList.size() >= 1 && timeIndex >= reuseInformationStartingFromTimeIndex )
  {
    assert( residualsTimeList.size() == solsTimeList.size() );
    residualsTimeList.push_front( residualsList );
    solsTimeList.push_front( solsList );
  }

  // Remove the last items from the residual list and solutions list
  // in order to ensure that at maximum nbReuse time steps
  // are included.
  while ( static_cast<int>( residualsTimeList.size() ) > nbReuse )
  {
    residualsTimeList.pop_back();
    solsTimeList.pop_back();
  }

  residualsList.clear();
  solsList.clear();

  assert( residualsTimeList.size() == solsTimeList.size() );
  assert( static_cast<int>( residualsTimeList.size() ) <= nbReuse );

  timeIndex++;
}

bool PostProcessing::isConvergence(
  const vector & xk,
  const vector & xkprev,
  bool residualCriterium
  )
{
  if ( residualCriterium )
    return fsi->isConvergence();

  return fsi->isConvergence( xk, xkprev );
}

void PostProcessing::iterationsConverged()
{
  // Save input/output information for next solve
  if ( residuals.size() >= 2 )
  {
    assert( residuals.size() == sols.size() );

    residualsList.push_front( residuals );
    solsList.push_front( sols );
  }

  residuals.clear();
  sols.clear();

  assert( residualsList.size() == solsList.size() );
}
