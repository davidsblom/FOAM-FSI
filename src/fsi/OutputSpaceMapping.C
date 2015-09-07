
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "OutputSpaceMapping.H"

using namespace fsi;

OutputSpaceMapping::OutputSpaceMapping(
    shared_ptr<SurrogateModel> fineModel,
    shared_ptr<SurrogateModel> surrogateModel,
    int maxIter,
    int maxUsedIterations,
    int nbReuse,
    int reuseInformationStartingFromTimeIndex,
    scalar singularityLimit,
    int order
    )
    :
    SpaceMapping( fineModel, surrogateModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ),
    surrogateModel( surrogateModel ),
    sols(),
    solsList(),
    solsTimeList(),
    order( order )
{
    assert( surrogateModel );
    assert( order == 0 || order == 1 || order == 2 );
}

OutputSpaceMapping::~OutputSpaceMapping()
{}

void OutputSpaceMapping::finalizeTimeStep()
{
    SpaceMapping::finalizeTimeStep();

    timeIndex--;

    // Save input/output information for next time step
    if ( nbReuse > 0 && solsList.size() >= 1 && timeIndex >= reuseInformationStartingFromTimeIndex )
        solsTimeList.push_front( solsList );

    timeIndex++;

    // Remove the last items from the residual list and solutions list
    // in order to ensure that at maximum nbReuse time steps
    // are included.
    while ( static_cast<int>( solsTimeList.size() ) > nbReuse )
        solsTimeList.pop_back();

    solsList.clear();

    assert( solsTimeList.size() == coarseResidualsTimeList.size() );
    assert( solsTimeList.size() == fineResidualsTimeList.size() );
    assert( static_cast<int>( solsTimeList.size() ) <= nbReuse );
}

void OutputSpaceMapping::iterationsConverged()
{
    SpaceMapping::iterationsConverged();

    // Save input/output information for next solve
    if ( sols.size() >= 2 )
        solsList.push_front( sols );

    sols.clear();

    assert( solsList.size() == coarseResidualsList.size() );
    assert( solsList.size() == fineResidualsList.size() );
}

void OutputSpaceMapping::performPostProcessing(
    const vector & y,
    const vector & x0,
    vector & xk,
    bool residualCriterium
    )
{
    assert( x0.rows() == xk.rows() );

    // Initialize variables

    int m = y.rows();
    int n = x0.rows();
    vector yk( m ), output( m ), R( m );
    xk = x0;
    vector xkprev = x0;
    output.setZero();
    R.setZero();
    coarseResiduals.clear();
    fineResiduals.clear();
    sols.clear();

    if ( timeIndex == 0 )
    {
        // Determine optimum of coarse model xstar
        if ( residualCriterium )
        {
            assert( y.norm() < 1.0e-14 );
            surrogateModel->optimize( x0, xk );
        }

        if ( !residualCriterium )
            surrogateModel->optimize( y, x0, xk );

        if ( !surrogateModel->allConverged() )
            Warning << "Output space mapping: surrogate model optimization process is not converged." << endl;
    }

    assert( xk.rows() == n );
    assert( x0.rows() == n );

    surrogateModel->setUseJacobian( false );

    // Initialize coarse model and fine model responses

    // Coarse model evaluation

    surrogateModel->evaluate( xk, output, R );
    coarseResiduals.push_back( R );
    assert( xk.rows() == n );
    assert( output.rows() == m );
    assert( R.rows() == m );

    // Fine model evaluation

    fineModel->evaluate( xk, output, R );
    fineResiduals.push_back( R );
    assert( output.rows() == m );
    assert( R.rows() == m );

    sols.push_back( xk );

    // Check convergence criteria
    if ( isConvergence( xk, xkprev, residualCriterium ) )
    {
        assert( fineModel->allConverged() );
        iterationsConverged();
        return;
    }

    assert( xk.rows() == n );

    for ( int k = 0; k < maxIter - 1; k++ )
    {
        // Determine the number of columns used to calculate the mapping matrix

        int nbCols = k;
        int nbColsCurrentTimeStep = nbCols;

        // Include information from previous optimization cycles

        for ( unsigned i = 0; i < fineResidualsList.size(); i++ )
            nbCols += fineResidualsList.at( i ).size() - 1;

        // Include information from previous time steps
        for ( unsigned i = 0; i < fineResidualsTimeList.size(); i++ )
            for ( unsigned j = 0; j < fineResidualsTimeList.at( i ).size(); j++ )
                nbCols += fineResidualsTimeList.at( i ).at( j ).size() - 1;

        // Update the design specification yk
        yk = coarseResiduals.at( k ) - (fineResiduals.at( k ) - y);

        if ( nbCols > 0 && order == 1 )
        {
            assert( fineResiduals.size() == coarseResiduals.size() );
            assert( sols.size() == fineResiduals.size() );
            assert( solsList.size() == fineResidualsList.size() );
            assert( solsList.size() == coarseResidualsList.size() );

            // Initialize mapping matrix
            matrix J = -fsi::matrix::Identity( m, m );
            fsi::vector d, dprev, deltad, deltax;
            int colIndex = 0;

            Info << "OSM(1) with ";
            Info << nbCols;
            Info << " cols for the Jacobian" << endl;

            // Include information from previous time steps

            for ( unsigned i = solsTimeList.size(); i-- > 0; )
            {
                for ( unsigned j = solsTimeList.at( i ).size(); j-- > 0; )
                {
                    for ( unsigned k = 0; k < solsTimeList.at( i ).at( j ).size() - 1; k++ )
                    {
                        colIndex++;

                        d = fineResidualsTimeList.at( i ).at( j ).at( k + 1 ) - coarseResidualsTimeList.at( i ).at( j ).at( k + 1 );
                        dprev = fineResidualsTimeList.at( i ).at( j ).at( k ) - coarseResidualsTimeList.at( i ).at( j ).at( k );
                        deltad = d - dprev;
                        deltax = solsTimeList.at( i ).at( j ).at( k + 1 ) - solsTimeList.at( i ).at( j ).at( k );

                        // Broyden update for the Jacobian matrix

                        if ( deltax.norm() < singularityLimit )
                            continue;

                        J += (deltad - J * deltax) / deltax.squaredNorm() * deltax.transpose();
                    }
                }
            }

            // Include information from previous optimization cycles

            for ( unsigned i = solsList.size(); i-- > 0; )
            {
                for ( unsigned j = 0; j < solsList.at( i ).size() - 1; j++ )
                {
                    colIndex++;

                    d = fineResidualsList.at( i ).at( j + 1 ) - coarseResidualsList.at( i ).at( j + 1 );
                    dprev = fineResidualsList.at( i ).at( j ) - coarseResidualsList.at( i ).at( j );
                    deltad = d - dprev;
                    deltax = solsList.at( i ).at( j + 1 ) - solsList.at( i ).at( j );

                    // Broyden update for the Jacobian matrix

                    if ( deltax.norm() < singularityLimit )
                        continue;

                    J += (deltad - J * deltax) / deltax.squaredNorm() * deltax.transpose();
                }
            }

            for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
            {
                colIndex++;

                d = fineResiduals.at( i + 1 ) - coarseResiduals.at( i + 1 );
                dprev = fineResiduals.at( i ) - coarseResiduals.at( i );
                deltad = d - dprev;
                deltax = sols.at( i + 1 ) - sols.at( i );

                // Broyden update for the Jacobian matrix

                if ( deltax.norm() < singularityLimit )
                    continue;

                J += (deltad - J * deltax) / deltax.squaredNorm() * deltax.transpose();
            }

            assert( colIndex == nbCols );

            surrogateModel->setUseJacobian( true );
            surrogateModel->setSurrogateData( xk, J );
        }

        if ( nbCols > 0 && order == 2 )
        {
            nbCols = std::min( nbCols, n );
            nbCols = std::min( nbCols, maxUsedIterations );
            nbColsCurrentTimeStep = std::min( nbCols, nbColsCurrentTimeStep );

            Info << "OSM(2) with ";
            Info << nbCols;
            Info << " cols for the Jacobian" << endl;

            matrix DeltaF( m, nbCols ), DeltaX( m, nbCols );

            int colIndex = 0;

            // Include information from current time step

            for ( int i = 0; i < nbColsCurrentTimeStep; i++ )
            {
                if ( colIndex >= DeltaF.cols() )
                    continue;

                DeltaF.col( colIndex ) = fineResiduals.back() - coarseResiduals.back();
                DeltaF.col( colIndex ) -= fineResiduals.at( k - 1 - i ) - coarseResiduals.at( k - 1 - i );
                DeltaX.col( colIndex ) = sols.back() - sols.at( k - 1 - i );
                colIndex++;
            }

            // Include information from previous optimization cycles

            for ( unsigned i = 0; i < solsList.size(); i++ )
            {
                for ( unsigned j = 0; j < solsList.at( i ).size() - 1; j++ )
                {
                    if ( colIndex >= DeltaF.cols() )
                        continue;

                    DeltaF.col( colIndex ) = fineResidualsList.at( i ).back() - coarseResidualsList.at( i ).back();

                    DeltaF.col( colIndex ) -= fineResidualsList.at( i ).at( fineResidualsList.at( i ).size() - 2 - j );
                    DeltaF.col( colIndex ) -= coarseResidualsList.at( i ).at( coarseResidualsList.at( i ).size() - 2 - j );

                    DeltaX.col( colIndex ) = solsList.at( i ).back() - solsList.at( i ).at( solsList.at( i ).size() - 2 - j );

                    colIndex++;
                }
            }

            // Include information from previous time steps
            for ( unsigned i = 0; i < solsTimeList.size(); i++ )
            {
                for ( unsigned j = 0; j < solsTimeList.at( i ).size(); j++ )
                {
                    for ( unsigned k = 0; k < solsTimeList.at( i ).at( j ).size() - 1; k++ )
                    {
                        assert( fineResidualsTimeList.at( i ).at( j ).size() == solsTimeList.at( i ).at( j ).size() );
                        assert( coarseResidualsTimeList.at( i ).at( j ).size() == solsTimeList.at( i ).at( j ).size() );

                        if ( colIndex >= DeltaF.cols() )
                            continue;

                        DeltaF.col( colIndex ) = fineResidualsTimeList.at( i ).at( j ).back();
                        DeltaF.col( colIndex ) -= coarseResidualsTimeList.at( i ).at( j ).back();

                        DeltaF.col( colIndex ) -= fineResidualsTimeList.at( i ).at( j ).at( fineResidualsTimeList.at( i ).at( j ).size() - 2 - k );
                        DeltaF.col( colIndex ) -= coarseResidualsTimeList.at( i ).at( j ).at( coarseResidualsTimeList.at( i ).at( j ).size() - 2 - k );

                        DeltaX.col( colIndex ) = solsTimeList.at( i ).at( j ).back();
                        DeltaX.col( colIndex ) -= solsTimeList.at( i ).at( j ).at( solsTimeList.at( i ).at( j ).size() - 2 - k );

                        colIndex++;
                    }
                }
            }

            assert( colIndex == nbCols );

            // Truncated singular value decomposition to solve for the
            // coefficients

            Eigen::JacobiSVD<matrix> svd( DeltaX, Eigen::ComputeThinU | Eigen::ComputeThinV );

            vector singularValues_inv = svd.singularValues();

            for ( unsigned int i = 0; i < singularValues_inv.rows(); ++i )
            {
                if ( svd.singularValues() ( i ) > singularityLimit )
                    singularValues_inv( i ) = 1.0 / svd.singularValues() ( i );
                else
                    singularValues_inv( i ) = 0;
            }

            matrix I = fsi::matrix::Identity( m, m );

            matrix pseudoDeltaX = svd.matrixV() * singularValues_inv.asDiagonal() * svd.matrixU().transpose();

            matrix J = -I + (DeltaF + DeltaX) * pseudoDeltaX;

            surrogateModel->setUseJacobian( true );
            surrogateModel->setSurrogateData( xk, J );
        }

        xkprev = xk;

        // Update the fine model optimum
        output.resize( n );
        surrogateModel->optimize( yk, xk, output );
        assert( output.rows() == n );
        assert( yk.rows() == m );
        assert( xk.rows() == n );

        if ( !surrogateModel->allConverged() )
            Warning << "Output space mapping: surrogate model optimization process is not converged." << endl;

        xk = output;

        surrogateModel->setUseJacobian( false );

        // Coarse model evaluation
        output.resize( m );
        surrogateModel->evaluate( xk, output, R );
        assert( xk.rows() == n );
        assert( output.rows() == m );
        assert( R.rows() == m );
        coarseResiduals.push_back( R );

        // Fine model evaluation

        fineModel->evaluate( xk, output, R );
        assert( xk.rows() == n );
        assert( output.rows() == m );
        assert( R.rows() == m );
        fineResiduals.push_back( R );

        sols.push_back( xk );

        // Check convergence criteria
        if ( isConvergence( xk, xkprev, residualCriterium ) )
        {
            assert( fineModel->allConverged() );
            iterationsConverged();
            break;
        }
    }
}

void OutputSpaceMapping::removeColumnFromMatrix(
    matrix & A,
    int col
    )
{
    for ( int j = col; j < A.cols() - 1; j++ )
        A.col( j ) = A.col( j + 1 );

    A.conservativeResize( A.rows(), A.cols() - 1 );
}
