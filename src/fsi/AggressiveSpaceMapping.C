
/*
 * Copyright [2016] <David Blom>
 */

#include "AggressiveSpaceMapping.H"

namespace fsi {

AggressiveSpaceMapping::AggressiveSpaceMapping(
    shared_ptr<SurrogateModel> fineModel,
    shared_ptr<SurrogateModel> coarseModel,
    int maxIter,
    int maxUsedIterations,
    int nbReuse,
    int reuseInformationStartingFromTimeIndex,
    scalar singularityLimit
    )
    :
    SpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit )
{}

AggressiveSpaceMapping::~AggressiveSpaceMapping()
{}

void AggressiveSpaceMapping::performPostProcessing(
    const vector & y,
    const vector & x0,
    vector & xk,
    bool residualCriterium
    )
{
    assert( x0.rows() == xk.rows() );

    // Initialize variables

    int m = y.rows();
    vector output( m ), R( m ), zstar( m ), zk( m ), xkprev( m );
    xk = x0;
    xkprev = x0;
    coarseResiduals.clear();
    fineResiduals.clear();

    // Determine optimum of coarse model zstar
    if ( residualCriterium )
    {
        assert( y.norm() < 1.0e-14 );
        coarseModel->optimize( x0, zstar );
    }

    if ( !residualCriterium )
        coarseModel->optimize( y, x0, zstar );

    if ( !coarseModel->allConverged() )
        Warning << "Surrogate model optimization process is not converged." << endl;

    xk = zstar;

    // Fine model evaluation

    fineModel->evaluate( xk, output, R );

    // Check convergence criteria
    if ( isConvergence( output, xk + y, residualCriterium ) )
    {
        iterationsConverged();
        return;
    }

    // Parameter extraction

    coarseModel->optimize( R - y, zstar, zk );

    if ( !coarseModel->allConverged() )
        Warning << "Surrogate model optimization process is not converged." << endl;

    coarseResiduals.push_back( zk );
    fineResiduals.push_back( xk );

    for ( int k = 0; k < maxIter - 1; k++ )
    {
        matrix J = fsi::matrix::Identity( m, m );

        // Determine the number of columns used to calculate the mapping matrix J

        int nbCols = k;
        int nbColsCurrentTimeStep = nbCols;

        // Include information from previous optimization cycles

        for ( auto && fineResiduals : fineResidualsList )
            nbCols += fineResiduals.size() - 1;

        // Include information from previous time steps

        for ( auto && fineResidualsList : fineResidualsTimeList )
            for ( auto && fineResiduals : fineResidualsList )
                nbCols += fineResiduals.size() - 1;

        if ( nbCols > 0 )
        {
            Info << "Aggressive space mapping with ";
            Info << nbCols;
            Info << " cols for the Jacobian" << endl;

            int colIndex = 0;

            // Include information from previous time steps

            for ( unsigned i = fineResidualsTimeList.size(); i-- > 0; )
            {
                for ( unsigned j = fineResidualsTimeList.at( i ).size(); j-- > 0; )
                {
                    assert( fineResidualsTimeList.at( i ).at( j ).size() >= 2 );
                    assert( coarseResidualsTimeList.at( i ).at( j ).size() >= 2 );

                    for ( unsigned k = 0; k < fineResidualsTimeList.at( i ).at( j ).size() - 1; k++ )
                    {
                        colIndex++;

                        vector deltax, deltaz;

                        deltax = fineResidualsTimeList.at( i ).at( j ).at( k + 1 ) - fineResidualsTimeList.at( i ).at( j ).at( k );
                        deltaz = coarseResidualsTimeList.at( i ).at( j ).at( k + 1 ) - coarseResidualsTimeList.at( i ).at( j ).at( k );

                        if ( deltax.norm() >= singularityLimit )
                        {
                            // Sherman–Morrison formula
                            J += (deltax - J * deltaz) / (deltax.transpose() * J * deltaz) * (deltax.transpose() * J);
                        }
                    }
                }
            }

            // Include information from previous optimization cycles

            for ( unsigned i = fineResidualsList.size(); i-- > 0; )
            {
                assert( fineResidualsList.at( i ).size() >= 2 );
                assert( coarseResidualsList.at( i ).size() >= 2 );

                for ( unsigned j = 0; j < fineResidualsList.at( i ).size() - 1; j++ )
                {
                    colIndex++;

                    vector deltax, deltaz;

                    deltax = fineResidualsList.at( i ).at( j + 1 ) - fineResidualsList.at( i ).at( j );
                    deltaz = coarseResidualsList.at( i ).at( j + 1 ) - coarseResidualsList.at( i ).at( j );

                    if ( deltax.norm() >= singularityLimit )
                    {
                        // Sherman–Morrison formula
                        J += (deltax - J * deltaz) / (deltax.transpose() * J * deltaz) * (deltax.transpose() * J);
                    }
                }
            }

            // Include information from previous iterations

            for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
            {
                colIndex++;

                vector deltax, deltaz;

                deltax = fineResiduals.at( i + 1 ) - fineResiduals.at( i );
                deltaz = coarseResiduals.at( i + 1 ) - coarseResiduals.at( i );

                if ( deltax.norm() >= singularityLimit )
                {
                    // Sherman–Morrison formula
                    J += (deltax - J * deltaz) / (deltax.transpose() * J * deltaz) * (deltax.transpose() * J);
                }
            }

            assert( colIndex == nbCols );
        }

        xkprev = xk;

        xk += J * (zstar - zk);

        // Fine model evaluation

        fineModel->evaluate( xk, output, R );

        // Check convergence criteria
        if ( isConvergence( output, xk + y, residualCriterium ) )
        {
            iterationsConverged();
            break;
        }

        // Parameter extraction

        coarseModel->optimize( R - y, zstar, zk );

        if ( !coarseModel->allConverged() )
            Warning << "Surrogate model optimization process is not converged." << endl;

        coarseResiduals.push_back( zk );
        fineResiduals.push_back( xk );
    }
}

} // namespace fsi
