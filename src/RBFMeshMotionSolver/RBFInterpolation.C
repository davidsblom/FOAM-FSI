
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolation.H"

namespace rbf
{
    Foam::debug::debugSwitch RBFInterpolation::debug( "RBFInterpolation", 0 );

    RBFInterpolation::RBFInterpolation( std::shared_ptr<RBFFunctionInterface> rbfFunction )
        :
        rbfFunction( rbfFunction ),
        polynomialTerm( true ),
        cpu( false ),
        computed( false ),
        n_A( 0 ),
        n_B( 0 ),
        dimGrid( 0 ),
        Hhat(),
        Phi(),
        lu(),
        fullPivLu(),
        llt(),
        positions(),
        positionsInterpolation()
    {
        assert( rbfFunction );
    }

    RBFInterpolation::RBFInterpolation(
        std::shared_ptr<RBFFunctionInterface> rbfFunction,
        bool polynomialTerm,
        bool cpu
        )
        :
        rbfFunction( rbfFunction ),
        polynomialTerm( polynomialTerm ),
        cpu( cpu ),
        computed( false ),
        n_A( 0 ),
        n_B( 0 ),
        dimGrid( 0 ),
        Hhat(),
        Phi(),
        lu(),
        fullPivLu(),
        llt(),
        positions(),
        positionsInterpolation()
    {
        assert( rbfFunction );
    }

    void RBFInterpolation::evaluateH(
        const matrix & positions,
        matrix & H
        )
    {
        // RBF function evaluation

        scalar r = 0;

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

        scalar r = 0;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_B; j++ )
            {
                r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
                Phi( j, i ) = rbfFunction->evaluate( r );
            }
        }
    }

    matrix RBFInterpolation::evaluatePhi(
        const matrix & positions,
        const vector & positionsInterpolation
        )
    {
        // Evaluate Phi which contains the evaluation of the radial basis function
        assert( positionsInterpolation.rows() == dimGrid );

        double r = 0;

        matrix rowPhi( 1, positions.rows() );

        for ( int i = 0; i < n_A; i++ )
        {
            r = ( positions.row( i ) - positionsInterpolation.transpose() ).norm();
            rowPhi( 0, i ) = rbfFunction->evaluate( r );
        }

        return rowPhi;
    }

    void RBFInterpolation::compute(
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        std::clock_t t = std::clock();

        if ( debug )
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
        matrix H( n_A, n_A ), Phi( n_B, n_A );

        if ( polynomialTerm )
        {
            H.resize( n_A + dimGrid + 1, n_A + dimGrid + 1 );
        }

        // Evaluate radial basis functions for matrix H

        evaluateH( positions, H );

        // Include polynomial contributions
        if ( polynomialTerm )
        {
            for ( int i = 0; i < n_A; i++ )
                H( n_A, i ) = 1;

            H.bottomLeftCorner( dimGrid, n_A ) = positions.block( 0, 0, n_A, dimGrid ).transpose();

            for ( int i = 0; i < dimGrid + 1; i++ )
                for ( int j = 0; j < dimGrid + 1; j++ )
                    H( H.rows() - dimGrid - 1 + i, H.rows() - dimGrid - 1 + j ) = 0;
        }

        if ( debug )
        {
            t = std::clock() - t;
            scalar runTime = static_cast<scalar>(t) / CLOCKS_PER_SEC;
            Info << "RBFInterpolation::debug 1. evaluate H = " << runTime << " s" << endl;
            t = std::clock();
        }

        if ( cpu )
        {
            this->positions = positions;
            this->positionsInterpolation = positionsInterpolation;

            // If the thin plate spline radial basis function is used,
            // use the LU decomposition to solve for the coefficients.
            // In case a wendland function is used, use the more efficient LLT
            // algorithm.
            // The LLT algorithm cannot be used for the thin plate spline function,
            // since the diagonal of matrix is zero for this function.
            std::shared_ptr<TPSFunction> function;
            function = std::dynamic_pointer_cast<TPSFunction>( rbfFunction );

            // THIJS: tmp switch off other solves to prevent crashing
            // if ( polynomialTerm )
            // fullPivLu.compute( H.selfadjointView<Eigen::Lower>() );
            // else
            // if ( function )
            // lu.compute( H.selfadjointView<Eigen::Lower>() );
            // else
            // llt.compute( H.selfadjointView<Eigen::Lower>() );
            fullPivLu.compute( H.selfadjointView<Eigen::Lower>() );
        }

        if ( not cpu )
        {
            // Evaluate Phi which contains the evaluation of the radial basis function

            if ( polynomialTerm )
                Phi.resize( n_B, n_A + dimGrid + 1 );

            evaluatePhi( positions, positionsInterpolation, Phi );

            // Include polynomial contributions in matrix Phi
            if ( polynomialTerm )
            {
                for ( int i = 0; i < Phi.rows(); i++ )
                    Phi( i, n_A ) = 1;

                Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );
            }

            if ( debug )
            {
                t = std::clock() - t;
                scalar runTime = static_cast<scalar>(t) / CLOCKS_PER_SEC;
                Info << "RBFInterpolation::debug 2. evaluate Phi = " << runTime << " s" << endl;
                t = std::clock();
            }

            // Compute the LU decomposition of the matrix H
            Eigen::FullPivLU<matrix> lu( H.selfadjointView<Eigen::Lower>() );

            // Compute interpolation matrix

            Hhat.noalias() = Phi * lu.inverse();

            Hhat.conservativeResize( n_B, n_A );

            if ( debug )
            {
                t = std::clock() - t;
                scalar runTime = static_cast<scalar>(t) / CLOCKS_PER_SEC;
                Info << "RBFInterpolation::debug 3. compute Hhat = " << runTime << " s" << endl;
            }
        }

        computed = true;
    }

    void RBFInterpolation::interpolate(
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
        if ( cpu && not computed )
            compute( positions, positionsInterpolation );

        assert( computed );

        std::clock_t t = std::clock();

        if ( cpu )
        {
            // If the thin plate spline radial basis function is used,
            // use the LU decomposition to solve for the coefficients.
            // In case a wendland function is used, use the more efficient LLT
            // algorithm.
            // The LLT algorithm cannot be used for the thin plate spline function,
            // since the diagonal of matrix is zero for this function.
            std::shared_ptr<TPSFunction> function;
            function = std::dynamic_pointer_cast<TPSFunction>( rbfFunction );

            matrix B, valuesLU( n_A, values.cols() ), Phi( n_B, n_A );

            if ( polynomialTerm )
            {
                Phi.resize( n_B, n_A + dimGrid + 1 );
                valuesLU.resize( n_A + dimGrid + 1, values.cols() );
            }

            valuesLU.setZero();
            valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;

            // THIJS: tmp switched off other solvers to prevent crashing
            // if ( polynomialTerm )
            // B = fullPivLu.solve( valuesLU );
            // else
            // if ( function )
            // B = lu.solve( valuesLU );
            // else
            // B = llt.solve( valuesLU );
            B = fullPivLu.solve( valuesLU );

            // === Building complete matrix === //

            /*evaluatePhi( positions, positionsInterpolation, Phi );
             *
             * if ( polynomialTerm )
             * {
             *  // Include polynomial contributions in matrix Phi
             *
             *  for ( int i = 0; i < Phi.rows(); i++ )
             *      Phi( i, n_A ) = 1;
             *
             *  Phi.topRightCorner( n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );
             * }
             *
             * valuesInterpolation.noalias() = Phi * B;*/

            // === Evaluating row by row ==== //
            if ( valuesInterpolation.rows() != n_B )
            {
                valuesInterpolation = matrix( n_B, values.cols() );
            }

            if ( polynomialTerm )
            {
                matrix rowPhi( 1, n_A + 1 + dimGrid );

                for ( int i = 0; i < n_B; i++ )
                {
                    rowPhi.topLeftCorner( 1, n_A ) = evaluatePhi( positions, positionsInterpolation.row( i ) );
                    rowPhi( 0, n_A ) = 1.0;
                    rowPhi.topRightCorner( 1, dimGrid ) = positionsInterpolation.row( i );
                    valuesInterpolation.row( i ) = rowPhi * B;
                }
            }
            else
            {
                for ( int i = 0; i < n_B; i++ )
                {
                    matrix rowPhi = evaluatePhi( positions, positionsInterpolation.row( i ) );
                    valuesInterpolation.row( i ) = rowPhi * B;
                }
            }
        }

        if ( not cpu )
        {
            // double ttmp = 0;

            // if ( debug == 3 )
            // {
            //     valuesInterpolation.setZero();
            //
            //     boolList bla( Pstream::nProcs(), true );
            //     reduce( bla, sumOp<boolList>() );
            //
            //     ttmp = std::clock();
            //
            //     // Pout << "Debug Interpolation of CPU " << Pstream::myProcNo() << " valuesInterpolation: " << valuesInterpolation.rows() << ", " << valuesInterpolation.cols() << nl;
            // }

            valuesInterpolation = Hhat * values;

            // valuesInterpolation.noalias() = Hhat * values;

            // if ( debug == 3 )
            // {
            //     ttmp = std::clock() - ttmp;
            //
            //     // Pout << "Debug Interpolation of CPU " << Pstream::myProcNo() << " of [" << Hhat.rows() << ", " << Hhat.cols() << "]x[" << values.rows() << ", " << values.cols() << "]: " << static_cast<float>(ttmp) / CLOCKS_PER_SEC << " s" << endl;
            // }
        }

        if ( debug )
        {
            t = std::clock() - t;
            scalar runTime = static_cast<scalar>(t) / CLOCKS_PER_SEC;
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
        const matrix & positions,
        const matrix & positionsInterpolation,
        const matrix & values,
        matrix & valuesInterpolation
        )
    {
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
            // THIJS: initialize Phi if empty
            if ( Phi.cols() == 0 )
            {
                Phi.conservativeResize( n_B, dimGrid + 1 );
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

        matrix valuesLU( H.rows(), values.cols() );
        valuesLU.setZero();
        valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;

        // THIJS: tmp change in solver to be used. When taking lots of point with compact supported function it fails.

        /*if ( polynomialTerm )
         * {
         *  fullPivLu.compute( H.selfadjointView<Eigen::Lower>() );
         *  B = fullPivLu.solve( valuesLU );
         * }
         * else
         * if ( function )
         * {
         *  lu.compute( H.selfadjointView<Eigen::Lower>() );
         *  B = lu.solve( valuesLU );
         * }
         * else
         * {
         *  llt.compute( H.selfadjointView<Eigen::Lower>() );
         *  B = llt.solve( valuesLU );
         * }*/

        fullPivLu.compute( H.selfadjointView<Eigen::Lower>() );
        B = fullPivLu.solve( valuesLU );

        // Evaluate Phi_BA which contains the evaluation of the radial basis function
        // This method is only used by the greedy algorithm, and the matrix Phi
        // is therefore enlarged at every greedy step.

        buildPhi( positions, positionsInterpolation );

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
        const matrix & positions,
        const matrix & positionsInterpolation
        )
    {
        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        int phiColsOld = Phi.cols();

        if ( polynomialTerm )
            Phi.conservativeResize( n_B, n_A + dimGrid + 1 );
        else
            Phi.conservativeResize( n_B, n_A );

        int nNewPoints = Phi.cols() - phiColsOld;

        if ( nNewPoints == Phi.cols() )
            nNewPoints = n_A;

        scalar r = 0;

        for ( int i = 0; i < nNewPoints; i++ )
        {
            int index = Phi.cols() - (i + 1);

            if ( polynomialTerm )
                index = Phi.cols() - 1 - dimGrid - (i + 1);

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

        matrix valuesLU( values.rows(), values.cols() );

        // resize valuesLU if polynomial is used
        if ( polynomialTerm )
        {
            valuesLU.conservativeResize( values.rows() + values.cols() + 1, values.cols() );
        }

        valuesLU.setZero(); // initialize all values zero

        // Set correct part of valuesLU equal to values
        if ( polynomialTerm )
        {
            valuesLU.topLeftCorner( values.rows(), values.cols() ) = values;
        }
        else
        {
            valuesLU = values;
        }

        matrix B;

        // THIJS: tmp change in solver to be used. When taking lots of point with compact supported function it fails.

        /*if ( polynomialTerm )
         * {
         *  B = fullPivLu.solve( valuesLU );
         * }
         * else
         * if ( function )
         * {
         *  B = lu.solve( valuesLU );
         * }
         * else
         * {
         *  B = llt.solve( valuesLU );
         * }*/
        B = fullPivLu.solve( valuesLU );

        valuesInterpolation.noalias() = Phi * B;

        assert( valuesInterpolation.rows() == n_B );
        assert( values.cols() == valuesInterpolation.cols() );
    }
}
