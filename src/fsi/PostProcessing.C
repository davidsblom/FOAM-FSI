
/*
 * Copyright [2016] <David Blom>
 */

#include "PostProcessing.H"

namespace fsi {

PostProcessing::PostProcessing(
    shared_ptr<MultiLevelFsiSolver> fsi,
    scalar initialRelaxation,
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
    residualsStageList(),
    solsStageList(),
    residualsTimeList(),
    solsTimeList(),
    k( 0 ),
    stageIndex( 0 ),
    initStage_( false )
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

void PostProcessing::finalizeStage()
{
    assert( initStage_ );
    assert( residualsStageList.size() == k - 1 );
    assert( residualsList.size() == solsList.size() );
    assert( residualsStageList.size() == solsStageList.size() );

    // Save input/output information for next stage
    if ( residualsList.size() >= 1 && timeIndex >= reuseInformationStartingFromTimeIndex )
    {
        // Push all the residuals into the residualsStageList at stageIndex
        assert( residualsList.size() == solsList.size() );

        for ( unsigned int i = 0; i < residualsList.size(); i++ )
        {
            residualsStageList.at( stageIndex ).push_front( residualsList.at( i ) );
            solsStageList.at( stageIndex ).push_front( solsList.at( i ) );
        }
    }

    residualsList.clear();
    solsList.clear();

    assert( residualsStageList.size() == solsStageList.size() );
    assert( residualsStageList.at( stageIndex ).size() == solsStageList.at( stageIndex ).size() );

    stageIndex++;
    initStage_ = false;
}

void PostProcessing::finalizeTimeStep()
{
    // Save input/output information for next time step
    if ( nbReuse > 0 && residualsStageList.size() >= 1 && timeIndex >= reuseInformationStartingFromTimeIndex )
    {
        assert( residualsTimeList.size() == solsTimeList.size() );
        residualsTimeList.push_front( residualsStageList );
        solsTimeList.push_front( solsStageList );
    }

    // Remove the last items from the residual list and solutions list
    // in order to ensure that at maximum nbReuse time steps
    // are included.
    while ( static_cast<int>( residualsTimeList.size() ) > nbReuse )
    {
        residualsTimeList.pop_back();
        solsTimeList.pop_back();
    }

    for ( unsigned int i = 0; i < residualsStageList.size(); i++ )
    {
        residualsStageList.at( i ).clear();
        solsStageList.at( i ).clear();
    }

    assert( residualsTimeList.size() == solsTimeList.size() );
    assert( static_cast<int>( residualsTimeList.size() ) <= nbReuse );

    timeIndex++;
}

void PostProcessing::initStage( int stageIndex )
{
    assert( not initStage_ );

    this->stageIndex = stageIndex;

    initStage_ = true;
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

void PostProcessing::iterationsConverged( bool keepIterations )
{
    // Save input/output information for next solve
    if ( residuals.size() >= 2 && keepIterations )
    {
        assert( residuals.size() == sols.size() );

        residualsList.push_front( residuals );
        solsList.push_front( sols );
    }

    residuals.clear();
    sols.clear();

    assert( residualsList.size() == solsList.size() );
}

void PostProcessing::setNumberOfImplicitStages( int k )
{
    this->k = k + 1;

    deque<deque<vector> > residualsList;

    residualsStageList.clear();
    solsStageList.clear();

    for ( int i = 0; i < k; i++ )
    {
        residualsStageList.push_front( residualsList );
        solsStageList.push_front( residualsList );
    }

    assert( static_cast<int>( residualsStageList.size() ) == k );
}

} // namespace fsi
