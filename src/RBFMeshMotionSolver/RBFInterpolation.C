
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"

namespace rbf
{
    RBFInterpolation::RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction )
        :
        rbfFunction( rbfFunction ),
        computed( false ),
        n_A( 0 ),
        n_B( 0 ),
        dimGrid( 0 ),
        Hhat()
    {
        assert( rbfFunction );
    }

    void RBFInterpolation::evaluateH(
        const matrix & positions,
        matrix & H
        )
    {
        // RBF function evaluation

        double r = 0;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = i; j < n_A; j++ )
            {
                r = ( positions.row( i ) - positions.row( j ) ).norm();
                H( j, i ) = rbfFunction->evaluate( r );
            }
        }
    }

    void RBFInterpolation::evaluatePhi(
        const matrix & positions,
        const matrix & positionsInterpolation,
        matrix & Phi
        )
    {
        // Evaluate Phi which contains the evaluation of the radial basis function

        double r = 0;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_B; j++ )
            {
                r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
                Phi( j, i ) = rbfFunction->evaluate( r );
            }
        }
    }

    void RBFInterpolation::compute(
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        std::clock_t t = std::clock();

        int debugLevel = Foam::debug::debugSwitch( "RBFInterpolation", 0 );

        if ( debugLevel > 0 )
        {
            Info << "RBFInterpolation::debug positions = " << positions.rows() << " x " << positions.cols() << endl;
            Info << "RBFInterpolation::debug positionsInterpolation = " << positionsInterpolation.rows() << " x " << positionsInterpolation.cols() << endl;
        }

        // Verify input

        assert( positions.cols() == positionsInterpolation.cols() );
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( positionsInterpolation.rows() > 0 );

        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        dimGrid = positions.cols();

        // Radial basis function interpolation
        // Initialize matrices H and Phi

        matrix H( n_A + dimGrid + 1, n_A + dimGrid + 1 ), Phi( n_B, n_A + dimGrid + 1 );

        // Evaluate radial basis functions for matrix H

        evaluateH( positions, H );

        // Include polynomial contributions

        for ( int i = 0; i < n_A; i++ )
            H( n_A, i ) = 1;

        H.bottomLeftCorner( dimGrid, n_A ) = positions.block( 0, 0, n_A, dimGrid ).transpose();

        for ( int i = 0; i < dimGrid + 1; i++ )
            for ( int j = 0; j < dimGrid + 1; j++ )
                H( H.rows() - dimGrid - 1 + i, H.rows() - dimGrid - 1 + j ) = 0;

        if ( debugLevel > 0 )
        {
            t = std::clock() - t;
            double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
            Info << "RBFInterpolation::debug 1. evaluate H = " << runTime << " s" << endl;
            t = std::clock();
        }

        // Evaluate Phi which contains the evaluation of the radial basis function

        evaluatePhi( positions, positionsInterpolation, Phi );

        // Include polynomial contributions in matrix Phi

        for ( int i = 0; i < Phi.rows(); i++ )
            Phi( i, n_A ) = 1;

        Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );

        if ( debugLevel > 0 )
        {
            t = std::clock() - t;
            double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
            Info << "RBFInterpolation::debug 2. evaluate Phi = " << runTime << " s" << endl;
            t = std::clock();
        }

        // Compute the LU decomposition of the matrix H
        //Eigen::PartialPivLU<matrix> lu( H.selfadjointView<Eigen::Lower>() );
        Eigen::FullPivLU<matrix> lu( H.selfadjointView<Eigen::Lower>() );

        // Compute interpolation matrix

        Hhat.noalias() = Phi * lu.inverse();

        Hhat.conservativeResize( n_B, n_A );

        if ( debugLevel > 0 )
        {
            t = std::clock() - t;
            double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
            Info << "RBFInterpolation::debug 3. compute Hhat = " << runTime << " s" << endl;
        }

        computed = true;
    }

    void RBFInterpolation::interpolate(
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        assert( computed );

        std::clock_t t = std::clock();

        valuesInterpolation.noalias() = Hhat * values;

        if ( Foam::debug::debugSwitch( "RBFInterpolation", 0 ) > 0 )
        {
            t = std::clock() - t;
            double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
            Info << "RBFInterpolation::debug 4. interpolation = " << runTime << " s" << endl;
        }

        assert( valuesInterpolation.rows() == n_B );
        assert( values.cols() == valuesInterpolation.cols() );
    }

    /*
     * Compute interpolation matrix and directly interpolate the values.
     * The algorithms solves for the coefficients, and explicitly
     * uses the coefficients to interpolate the data to the new positions.
     */
    void RBFInterpolation::interpolate(
        bool polynomialTerm,
        const matrix & positions,
        const matrix & positionsInterpolation,
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        // Verify input

        assert( positions.cols() == positionsInterpolation.cols() );
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( positionsInterpolation.rows() > 0 );

        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        dimGrid = positions.cols();

        // Radial basis function interpolation

        // Initialize matrices

        matrix H( n_A, n_A );

        if ( polynomialTerm )
            H.resize( n_A + dimGrid + 1, n_A + dimGrid + 1 );

        // RBF function evaluation

        evaluateH( positions, H );

        // Include polynomial contributions

        if ( polynomialTerm )
        {
            //THIJS: initialize Phi if empty
            if(Phi.cols() == 0)
            {
                Phi.conservativeResize( n_B, dimGrid + 1);
            }

            for ( int i = 0; i < n_A; i++ )
                H( n_A, i ) = 1;

            H.bottomLeftCorner( dimGrid, n_A ) = positions.block( 0, 0, n_A, dimGrid ).transpose();

            for ( int i = 0; i < dimGrid + 1; i++ )
                for ( int j = 0; j < dimGrid + 1; j++ )
                    H( H.rows() - dimGrid - 1 + i, H.rows() - dimGrid - 1 + j ) = 0;
        }

        // Calculate coefficients gamma and beta

        // If the thin plate spline radial basis function is used,
        // use the LU decomposition to solve for the coefficients.
        // In case a wendland function is used, use the more efficient LLT
        // algorithm.
        // The LLT algorithm cannot be used for the thin plate spline function,
        // since the diagonal of matrix is zero for this function.
        std::shared_ptr<TPSFunction> function;
        function = std::dynamic_pointer_cast<TPSFunction>( rbfFunction );

        matrix B;

        matrix valuesLU ( H.rows(), values.cols() );
        valuesLU.setZero();
        valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;

        if ( polynomialTerm )
        {
            fullPivLu.compute( H.selfadjointView<Eigen::Lower>() );
            B = fullPivLu.solve( valuesLU );
        }
        else
        if ( function )
        {
            lu.compute( H.selfadjointView<Eigen::Lower>() );
            B = lu.solve( valuesLU );
        }
        else
        {
            llt.compute( H.selfadjointView<Eigen::Lower>() );
            B = llt.solve( valuesLU );
        }

        // Evaluate Phi_BA which contains the evaluation of the radial basis function
        // This method is only used by the greedy algorithm, and the matrix Phi
        // is therefore enlarged at every greedy step.

        buildPhi( polynomialTerm, positions, positionsInterpolation );

        if ( polynomialTerm )
        {
            // Include polynomial contributions in matrix Phi

            for ( int i = 0; i < Phi.rows(); i++ )
                Phi( i, n_A ) = 1;

            Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );
        }

        valuesInterpolation.noalias() = Phi * B;

        computed = true;
    }

    void RBFInterpolation::buildPhi(
        bool polynomialTerm,
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        int phiColsOld = Phi.cols();

        if ( polynomialTerm )
            Phi.conservativeResize(n_B, n_A + dimGrid + 1 );
        else
            Phi.conservativeResize( n_B, n_A );

        int nNewPoints = Phi.cols() - phiColsOld;

        double r = 0;
        for (int i = 0; i < nNewPoints; i++)
        {
            int index = Phi.cols() - (i+1);
            if ( polynomialTerm )
                index = Phi.cols() - 1 - dimGrid - (i+1);

            for ( int j = 0; j < n_B; j++ )
            {

                r = ( positions.row( index ) - positionsInterpolation.row( j ) ).norm();
                Phi( j, index ) = rbfFunction->evaluate( r );
            }
        }
    }

    /*
    * This function is only called by the RBFCoarsening class.
    * It is assumed that the polynomial term is included in the
    * interpolation, and that the fullPivLu decomposition is
    * used to solve for the coefficients B.
    */
    void RBFInterpolation::interpolate2(
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        assert( computed );

        // If the thin plate spline radial basis function is used,
        // use the LU decomposition to solve for the coefficients.
        // In case a wendland function is used, use the more efficient LLT
        // algorithm.
        // The LLT algorithm cannot be used for the thin plate spline function,
        // since the diagonal of matrix is zero for this function.
        std::shared_ptr<TPSFunction> function;
        function = std::dynamic_pointer_cast<TPSFunction>( rbfFunction );

        matrix valuesLU ( values.rows() + values.cols() + 1, values.cols() );
        valuesLU.setZero();
        valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;

        matrix B = fullPivLu.solve( valuesLU );

        valuesInterpolation.noalias() = Phi * B;

        assert( valuesInterpolation.rows() == n_B );
        assert( values.cols() == valuesInterpolation.cols() );
    }
}
