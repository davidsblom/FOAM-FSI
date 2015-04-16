
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
    int maxUsedIterations,
    int nbReuse,
    int reuseInformationStartingFromTimeIndex,
    double singularityLimit
    )
    :
    fineModel( fineModel ),
    coarseModel( coarseModel ),
    maxIter( maxIter ),
    maxUsedIterations( maxUsedIterations ),
    nbReuse( nbReuse ),
    reuseInformationStartingFromTimeIndex( reuseInformationStartingFromTimeIndex ),
    singularityLimit( singularityLimit ),
    timeIndex( 0 ),
    coarseResiduals(),
    fineResiduals(),
    coarseResidualsList(),
    fineResidualsList(),
    coarseResidualsTimeList(),
    fineResidualsTimeList()
{
    assert( fineModel );
    assert( coarseModel );
    assert( maxIter > 0 );
    assert( nbReuse >= 0 );
    assert( reuseInformationStartingFromTimeIndex >= 0 );
    assert( singularityLimit > 0 );
    assert( singularityLimit < 1 );
    assert( maxUsedIterations > 0 );
};

void SpaceMapping::finalizeTimeStep()
{
    // Save input/output information for next time step
    if ( nbReuse > 0 && fineResidualsList.size() >= 1 && timeIndex >= reuseInformationStartingFromTimeIndex )
    {
        assert( fineResidualsTimeList.size() == coarseResidualsTimeList.size() );
        fineResidualsTimeList.push_front( fineResidualsList );
        coarseResidualsTimeList.push_front( coarseResidualsList );
    }

    // Remove the last items from the residual list and solutions list
    // in order to ensure that at maximum nbReuse time steps
    // are included.
    while ( static_cast<int>( fineResidualsTimeList.size() ) > nbReuse )
    {
        fineResidualsTimeList.pop_back();
        coarseResidualsTimeList.pop_back();
    }

    fineResidualsList.clear();
    coarseResidualsList.clear();

    assert( fineResidualsTimeList.size() == coarseResidualsTimeList.size() );
    assert( coarseResiduals.size() == fineResiduals.size() );
    assert( static_cast<int>( fineResidualsTimeList.size() ) <= nbReuse );

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
    // Save input/output information for next solve
    if ( fineResiduals.size() >= 2 )
    {
        assert( fineResiduals.size() == coarseResiduals.size() );

        fineResidualsList.push_front( fineResiduals );
        coarseResidualsList.push_front( coarseResiduals );
    }

    fineResiduals.clear();
    coarseResiduals.clear();

    assert( fineResidualsList.size() == coarseResidualsList.size() );
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
