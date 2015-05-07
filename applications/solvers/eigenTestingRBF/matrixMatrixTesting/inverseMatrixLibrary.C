
/*
* Author
*   David Blom, TU Delft. All rights reserved.
*/

#include "inverseMatrixLibrary.H"
#include <ctime>
#include <iostream>

namespace rbf
{
    // ===================== Constructor ================== //
    inverseMatrixLibrary::inverseMatrixLibrary( std::shared_ptr<RBFFunctionInterface> rbfFunction )
    :
    rbfFunction( rbfFunction ),
    n_A( 0 ),
    dimGrid( 0 ),
    C()
    {
        assert( rbfFunction );
    }

    void inverseMatrixLibrary::updateInverse(
        const matrix & oldCinverse,
        const matrix & newpoint,
        const matrix & oldPositions,
        matrix & Cinverse
    )
    {
        clock_t t0 = std::clock();

        // Verify input
        assert( newpoint.rows()==1);
        assert( oldPositions.rows() > 0 );
        assert( oldPositions.cols() > 0 );
        assert( newpoint.cols()==oldPositions.cols());
        assert( oldCinverse.cols() == oldCinverse.rows() );
        assert( oldCinverse.rows() == oldPositions.rows());
        assert( Cinverse.cols() == Cinverse.rows() );
        assert( Cinverse.rows() == oldCinverse.rows() + 1);

        n_A = oldPositions.rows();
        //dimGrid = oldPositions.cols();

        // Initialize variables
        matrix newCol(n_A, 1);
        newCol.setZero();

        // RBF function evaluation
        double r;
        for ( int i = 0; i < n_A; i++ )
        {
            r = (oldPositions.row(i) - newpoint.row(0)).norm();
            newCol(i,0) = rbfFunction->evaluate(r);
        }
        //std::cout << "Created new col/row" << std::endl;

        clock_t t1 = std::clock();

        //Construct elements for update
        matrix D(1,1);
        D(0,0) = rbfFunction->evaluate(0.0);
        matrix E = oldCinverse*newCol;
        matrix Xinv = (D - newCol.transpose()*E);
        matrix X(1,1);
        X(0,0) = 1.0/Xinv(0,0);
        //std::cout << "Created all variables needed for update" << std::endl;

        //Make update
        Cinverse.topLeftCorner(n_A,n_A) = oldCinverse + X(0,0)*(E*E.transpose());
        Cinverse.topRightCorner(n_A,1) = -X(0,0)*E;
        Cinverse.bottomLeftCorner(1,n_A) = -X(0,0)*E.transpose();
        Cinverse.bottomRightCorner(1,1) = X;
        clock_t t2 = std::clock();

        //Timings
        double t01 = 1000.0 * (t1-t0)/CLOCKS_PER_SEC;
        double t12 = 1000.0 * (t2-t1)/CLOCKS_PER_SEC;
        std::cout << "Timings: build new col|update inverse = " << t01 << "|" << t12 << std::endl;
    }

    // ========================================================== //
    // = Compute inverse using eigen
    // ========================================================== //
    void inverseMatrixLibrary::computeInverseEigen(
        const matrix & positions,
        matrix & Cinverse
    )
    {
        clock_t t0 = std::clock();
        // Verify input
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( Cinverse.cols() == Cinverse.rows() );
        assert( Cinverse.rows() == positions.rows());

        n_A = positions.rows();
        dimGrid = positions.cols();

        // Radial basis function interpolation

        // Initialize variables
        matrix Q_A( n_A, dimGrid + 1 );
        Q_A.setZero();

        // Initialize sparse matrices
        matrix C( n_A + dimGrid + 1, n_A + dimGrid + 1 );
        C.setZero();

        // RBF function evaluation
        double r;
        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_A; j++ )
            {
                r = ( positions.row( i ) - positions.row( j ) ).norm();
                C( j, i ) = rbfFunction->evaluate( r );
            }
        }

        // Evaluate Q_A
        for ( int i = 0; i < n_A; i++ )
        Q_A( i, 0 ) = 1;

        Q_A.block( 0, 1, n_A, dimGrid ) = positions.block( 0, 0, n_A, dimGrid );

        // Build the matrix H

        C.topRightCorner( Q_A.rows(), Q_A.cols() ) = Q_A;
        C.bottomLeftCorner( Q_A.cols(), Q_A.rows() ) = Q_A.transpose();

        // Compute the LU decomposition of the matrix H
        clock_t t1 = std::clock();
        Eigen::FullPivLU<matrix> lu( C );
        clock_t t2 = std::clock();
        Cinverse = lu.inverse();
        clock_t t3 = std::clock();
        //matrix betas = C.fullPivLu().solve( values );

        double t01 = 1000.0 * (t1-t0)/CLOCKS_PER_SEC;
        double t12 = 1000.0 * (t2-t1)/CLOCKS_PER_SEC;
        double t23 = 1000.0 * (t3-t2)/CLOCKS_PER_SEC;
        //std::cout << "Timings: build C|LU|inverse = " << t01 << "|" << t12 << "|" << t23 << std::endl;
    }

    // ========================================================== //
    // = Compute inverse using eigen
    // ========================================================== //
    void inverseMatrixLibrary::solveSystemEigen(
        const matrix & positions,
        const matrix & displacement,
        matrix & solution
    )
    {
        clock_t t0 = std::clock();
        // Verify input
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( displacement.cols() == 1);
        assert( displacement.rows() == positions.rows() );
        assert( solution.rows() == displacement.rows() );

        n_A = positions.rows();
        dimGrid = positions.cols();

        // Radial basis function interpolation

        // Initialize variables
        matrix Q_A( n_A, dimGrid + 1 );
        Q_A.setZero();

        // Initialize sparse matrices
        matrix C( n_A + dimGrid + 1, n_A + dimGrid + 1 );
        C.setZero();

        // RBF function evaluation
        double r;
        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_A; j++ )
            {
                r = ( positions.row( i ) - positions.row( j ) ).norm();
                C( j, i ) = rbfFunction->evaluate( r );
            }
        }

        // Evaluate Q_A
        for ( int i = 0; i < n_A; i++ )
        Q_A( i, 0 ) = 1;

        Q_A.block( 0, 1, n_A, dimGrid ) = positions.block( 0, 0, n_A, dimGrid );

        // Build the matrix C
        C.topRightCorner( Q_A.rows(), Q_A.cols() ) = Q_A;
        C.bottomLeftCorner( Q_A.cols(), Q_A.rows() ) = Q_A.transpose();

        //build rhs
        matrix rhs( n_A + dimGrid + 1, 1);
        rhs.setZero();
        rhs.topLeftCorner(n_A,1) = displacement;

        // Solve the system
        clock_t t1 = std::clock();
        solution = C.fullPivLu().solve( rhs );
        clock_t t2 = std::clock();

        double t01 = 1000.0 * (t1-t0)/CLOCKS_PER_SEC;
        double t12 = 1000.0 * (t2-t1)/CLOCKS_PER_SEC;
        //std::cout << "Timings: build C|solve = " << t01 << "|" << t12 << std::endl;
    }

    /*void inverseMatrixLibrary::compute(
        const matrix & positions,
        const matrix & positionsInterpolation
    )
    {
        // Verify input

        assert( positions.cols() == positionsInterpolation.cols() );
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( positionsInterpolation.rows() > 0 );

        n_A = positions.rows();
        n_B = positionsInterpolation.rows();
        Hhat.resize( n_B, n_A );
        Hhat.setZero();
        dimGrid = positions.cols();

        // Radial basis function interpolation

        // Initialize variables

        matrix Q_A( n_A, dimGrid + 1 ), Q_B( n_B, dimGrid + 1 );
        Q_A.setZero();
        Q_B.setZero();

        // Initialize sparse matrices

        matrix H( n_A + dimGrid + 1, n_A + dimGrid + 1 ), Phi( n_B, n_A + dimGrid + 1 ), Phi_AA( n_A, n_A );
        H.setZero();

        // RBF function evaluation

        double r;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_A; j++ )
            {
                r = ( positions.row( i ) - positions.row( j ) ).norm();
                H( j, i ) = rbfFunction->evaluate( r );
            }
        }

        // Evaluate Q_A
        for ( int i = 0; i < n_A; i++ )
        Q_A( i, 0 ) = 1;

        Q_A.block( 0, 1, n_A, dimGrid ) = positions.block( 0, 0, n_A, dimGrid );

        // Build the matrix H

        H.topRightCorner( Q_A.rows(), Q_A.cols() ) = Q_A;
        H.bottomLeftCorner( Q_A.cols(), Q_A.rows() ) = Q_A.transpose();

        // Evaluate Phi_BA which contains the evaluation of the radial basis function

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_B; j++ )
            {
                r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
                Phi( j, i ) = rbfFunction->evaluate( r );
            }
        }

        // Evaluate Q_B

        for ( int i = 0; i < n_B; i++ )
        Q_B( i, 0 ) = 1;

        Q_B.block( 0, 1, n_B, dimGrid ) = positionsInterpolation.block( 0, 0, n_B, dimGrid );

        // Determine the matrix Phi

        Phi.topRightCorner( Q_B.rows(), Q_B.cols() ) = Q_B;

        // Compute the LU decomposition of the matrix H

        Eigen::FullPivLU<matrix> lu( H );
        matrix Hinverse = lu.inverse();

        // Matrix matrix multiplication
        // Hhat = Phi * inv(H)

        Hhat.setZero();

        for ( int i = 0; i < Hhat.rows(); i++ )
        for ( int j = 0; j < Hhat.cols(); j++ )
        for ( int k = 0; k < Phi.cols(); k++ )
        Hhat( i, j ) += Phi( i, k ) * Hinverse( k, j );

        computed = true;
    }

    void inverseMatrixLibrary::interpolate(
        const matrix & values,
        matrix & valuesInterpolation
    )
    {
        assert( computed );

        valuesInterpolation = Hhat * values;

        assert( valuesInterpolation.rows() == n_B );
        assert( values.cols() == valuesInterpolation.cols() );
    }

    void inverseMatrixLibrary::interpolate(
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
        Hhat.resize( n_B, n_A );
        Hhat.setZero();
        dimGrid = positions.cols();

        // Radial basis function interpolation

        // Initialize matrices

        matrix H( n_A, n_A ), Phi_AA( n_A, n_A );
        H.setZero();
        Phi.conservativeResize( n_B, n_A );

        // RBF function evaluation

        double r;

        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_A; j++ )
            {
                r = ( positions.row( i ) - positions.row( j ) ).norm();
                H( j, i ) = rbfFunction->evaluate( r );
            }
        }

        // Calculate coefficients gamma and beta

        matrix B = H.fullPivLu().solve( values );

        // Evaluate Phi_BA which contains the evaluation of the radial basis function
        // This method is only used by the greedy algorithm, and the matrix Phi
        // is therefore enlarged at every greedy step.

        int i = Phi.cols() - 1;

        for ( int j = 0; j < n_B; j++ )
        {
            r = ( positions.row( i ) - positionsInterpolation.row( j ) ).norm();
            Phi( j, i ) = rbfFunction->evaluate( r );
        }

        valuesInterpolation = Phi * B;
    }*/

    void inverseMatrixLibrary::getControlMatrix(
        const matrix & positions,
        matrix & Crbf
    )
    {
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( Crbf.cols() == Crbf.rows() );
        assert( Crbf.rows() == positions.rows());

        // RBF function evaluation
        double r;
        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < n_A; j++ )
            {
                r = ( positions.row(i) - positions.row(j)).norm();
                Crbf(j,i) = rbfFunction->evaluate(r);
            }
        }
    }

    void inverseMatrixLibrary::getControlMatrixLower(
        const matrix & positions,
        matrix & Crbf
    )
    {
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( Crbf.cols() == Crbf.rows() );
        assert( Crbf.rows() == positions.rows());

        // RBF function evaluation
        double r;
        for ( int i = 0; i < n_A; i++ )
        {
            for ( int j = 0; j < i+1; j++ )
            {
                r = ( positions.row(i) - positions.row(j)).norm();
                Crbf(i,j) = rbfFunction->evaluate(r);
            }
        }
    }

    void inverseMatrixLibrary::getControlMatrixUpper(
        const matrix & positions,
        matrix & Crbf
    )
    {
        assert( positions.rows() > 0 );
        assert( positions.cols() > 0 );
        assert( Crbf.cols() == Crbf.rows() );
        assert( Crbf.rows() == positions.rows());

        // RBF function evaluation
        double r;
        for (int i=0; i<n_A; i++)
        {
            for (int j=i; j<n_A; j++ )
            {
                r = ( positions.row(i) - positions.row(j)).norm();
                Crbf(i,j) = rbfFunction->evaluate(r);
            }
        }
    }
}
