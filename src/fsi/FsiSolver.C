
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef EIGEN_CONFIG_H_
 #define EIGEN_CONFIG_H_

 #include <boost/serialization/array.hpp>
 #define EIGEN_DENSEBASE_PLUGIN "EigenDenseBaseAddons.h"

 #include <Eigen/Core>

#endif  // EIGEN_CONFIG_H_

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "FsiSolver.H"

/// Boost Serialization Helper

template <typename T>
bool serialize(
    const T & data,
    const std::string & filename
    )
{
    std::ofstream ofs( filename.c_str(), std::ios::out );

    if ( !ofs.is_open() )
        return false;

    {
        boost::archive::binary_oarchive oa( ofs );
        oa << data;
    }
    ofs.close();
    return true;
}

template <typename T>
bool deSerialize(
    T & data,
    const std::string & filename
    )
{
    std::ifstream ifs( filename.c_str(), std::ios::in );

    if ( !ifs.is_open() )
        return false;

    {
        boost::archive::binary_iarchive ia( ifs );
        ia >> data;
    }
    ifs.close();
    return true;
}

namespace fsi
{
    FsiSolver::FsiSolver(
        shared_ptr<BaseMultiLevelSolver> fluid,
        shared_ptr<BaseMultiLevelSolver> solid,
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures,
        bool parallel,
        int extrapolationOrder
        )
        :
        init( false ),
        fluid( fluid ),
        solid( solid ),
        convergenceMeasures( convergenceMeasures ),
        parallel( parallel ),
        N( fluid->N ),
        iter( 0 ),
        nbIter( 0 ),
        allConverged( false ),
        x(),
        extrapolationOrder( extrapolationOrder ),
        previousSolutions(),
        saveState( false ),
        saveInterval( 0 )
    {
        // Verify input parameters
        assert( fluid );
        assert( solid );
        assert( N > 0 );
        assert( convergenceMeasures );
        assert( convergenceMeasures->size() > 0 );
        assert( this->fluid->N == N );
        assert( fluid->N > 0 );
        assert( solid->N > 0 );
        assert( extrapolationOrder >= 0 );
        assert( extrapolationOrder <= 2 );

        // Initialize control variable x
        if ( parallel )
        {
            x.resize( solid->data.rows() * solid->data.cols() + fluid->data.rows() * fluid->data.cols() );

            for ( int i = 0; i < solid->data.cols(); i++ )
                x.segment( i * solid->data.rows(), solid->data.rows() ) = solid->data.col( i );

            for ( int j = 0; j < fluid->data.cols(); j++ )
                x.segment( j * fluid->data.rows() + solid->data.rows() * solid->data.cols(), fluid->data.rows() ) = fluid->data.col( j );
        }

        if ( !parallel )
        {
            x.resize( solid->data.rows() * solid->data.cols() );

            for ( int i = 0; i < solid->data.cols(); i++ )
                x.segment( i * solid->data.rows(), solid->data.rows() ) = solid->data.col( i );
        }

        // Save initial solution for extrapolation
        previousSolutions.push_front( x );
    }

    FsiSolver::FsiSolver(
        shared_ptr<BaseMultiLevelSolver> fluid,
        shared_ptr<BaseMultiLevelSolver> solid,
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures,
        bool parallel,
        int extrapolationOrder,
        bool saveState,
        int saveInterval
        )
        :
        init( false ),
        fluid( fluid ),
        solid( solid ),
        convergenceMeasures( convergenceMeasures ),
        parallel( parallel ),
        N( fluid->N ),
        iter( 0 ),
        nbIter( 0 ),
        allConverged( false ),
        x(),
        extrapolationOrder( extrapolationOrder ),
        previousSolutions(),
        saveState( saveState ),
        saveInterval( saveInterval )
    {
        // Verify input parameters
        assert( fluid );
        assert( solid );
        assert( N > 0 );
        assert( convergenceMeasures );
        assert( convergenceMeasures->size() > 0 );
        assert( this->fluid->N == N );
        assert( fluid->N > 0 );
        assert( solid->N > 0 );
        assert( extrapolationOrder >= 0 );
        assert( extrapolationOrder <= 2 );
        assert( saveInterval >= 0 );

        // Initialize control variable x
        if ( parallel )
        {
            x.resize( solid->data.rows() * solid->data.cols() + fluid->data.rows() * fluid->data.cols() );

            for ( int i = 0; i < solid->data.cols(); i++ )
                x.segment( i * solid->data.rows(), solid->data.rows() ) = solid->data.col( i );

            for ( int j = 0; j < fluid->data.cols(); j++ )
                x.segment( j * fluid->data.rows() + solid->data.rows() * solid->data.cols(), fluid->data.rows() ) = fluid->data.col( j );
        }

        if ( !parallel )
        {
            x.resize( solid->data.rows() * solid->data.cols() );

            for ( int i = 0; i < solid->data.cols(); i++ )
                x.segment( i * solid->data.rows(), solid->data.rows() ) = solid->data.col( i );
        }

        // Save initial solution for extrapolation
        previousSolutions.push_front( x );
    }

    FsiSolver::~FsiSolver()
    {}

    /*
     * Evaluate the FSI residual. Two residuals are defined,
     * following from the fixed point equation x = S o F (x),
     * and derived from the parallel fixed point equations:
     * x = S ( p ), p = F ( x ).
     */
    void FsiSolver::evaluate(
        const vector & input,
        vector & output,
        vector & R
        )
    {
        assert( fluid->init );
        assert( solid->init );
        assert( input.rows() == R.rows() );
        assert( input.rows() == output.rows() );
        assert( input.rows() == N || input.rows() == 2 * N );
        assert( input.rows() == x.rows() );

        Info << endl << "Time = " << fluid->t << ", iteration: " << iter + 1 << endl;

        matrix a = input.head( N );

        if ( !parallel )
        {
            assert( input.rows() == N );

            // Initialize variables
            matrix p( N, 1 ), aout( N, 1 );

            fluid->solve( a, p );
            solid->solve( p, aout );

            output = aout.col( 0 );

            fluid->couplingData.setDataOld( fluid->couplingData.data );
            fluid->couplingData.setData( p );
            solid->couplingData.setDataOld( a );
            solid->couplingData.setData( aout );
        }

        if ( parallel )
        {
            assert( input.rows() == 2 * N );

            // Initialize variables
            matrix p( N, 1 ), pout( N, 1 ), aout( N, 1 );

            p = input.tail( N );

            fluid->solve( a, pout );
            solid->solve( p, aout );

            output.head( N ) = aout.col( 0 );
            output.tail( N ) = pout.col( 0 );

            fluid->couplingData.setDataOld( p );
            fluid->couplingData.setData( pout );
            solid->couplingData.setDataOld( a );
            solid->couplingData.setData( aout );
        }

        // Calculate residual
        R = output - input;

        // Increment iterators
        iter++;
        nbIter++;
    }

    void FsiSolver::extrapolateData()
    {
        assert( !init );

        if ( extrapolationOrder == 0 )
            return;

        if ( previousSolutions.size() == 1 )
            return;

        bool firstOrderExtrapolation = previousSolutions.size() == 2 || (previousSolutions.size() == 3 && extrapolationOrder == 1);
        bool secondOrderExtrapolation = extrapolationOrder == 2 && previousSolutions.size() == 3;

        if ( firstOrderExtrapolation )
        {
            Info << "Perform first order extrapolation" << endl;

            assert( !secondOrderExtrapolation );
            assert( previousSolutions.size() == 2 || previousSolutions.size() == 3 );
            assert( previousSolutions.at( 0 ).rows() == previousSolutions.at( 1 ).rows() );

            x = 2.0 * previousSolutions.at( 0 ).array() - previousSolutions.at( 1 ).array();
        }
        else
        if ( secondOrderExtrapolation )
        {
            Info << "Perform second order extrapolation" << endl;

            assert( !firstOrderExtrapolation );
            assert( previousSolutions.size() == 3 );
            assert( previousSolutions.at( 0 ).rows() == previousSolutions.at( 1 ).rows() );
            assert( previousSolutions.at( 0 ).rows() == previousSolutions.at( 2 ).rows() );

            x = 2.5 * previousSolutions.at( 0 ).array() - 2.0 * previousSolutions.at( 1 ).array() + 0.5 * previousSolutions.at( 2 ).array();
        }

        assert( (firstOrderExtrapolation && !secondOrderExtrapolation) || (!firstOrderExtrapolation && secondOrderExtrapolation) );
    }

    void FsiSolver::finalizeTimeStep()
    {
        assert( init );

        fluid->couplingData.finalizeTimeStep();
        solid->couplingData.finalizeTimeStep();

        fluid->finalizeTimeStep();
        solid->finalizeTimeStep();

        // Save solution for extrapolation of solutions to next time step
        assert( previousSolutions.at( 0 ).rows() == x.rows() );
        previousSolutions.push_front( x );

        // Only save solutions of three time steps
        while ( previousSolutions.size() > 3 )
            previousSolutions.pop_back();

        assert( previousSolutions.size() <= 3 );

        // Write current state to disk
        // if ( saveState && saveInterval % fluid->timeIndex == 0 && Pstream::myProcNo() == 0 )
        if ( Pstream::myProcNo() == 0 )
        {
            serialize( x, "fsi-interface-" + std::to_string( fluid->timeIndex ) + ".bin" );
        }

        init = false;
    }

    void FsiSolver::initTimeStep()
    {
        assert( !init );
        assert( !fluid->init );
        assert( !solid->init );

        iter = 0;
        allConverged = false;

        fluid->initTimeStep();
        solid->initTimeStep();

        extrapolateData();

        init = true;
        assert( fluid->init );
        assert( solid->init );
        assert( fluid->t > 0 );
    }

    bool FsiSolver::isConvergence()
    {
        assert( false );
        return true;
    }

    bool FsiSolver::isConvergence(
        const vector &,
        const vector &
        )
    {
        assert( false );
        return true;
    }

    bool FsiSolver::isRunning()
    {
        assert( !init );

        return fluid->isRunning();
    }

    void FsiSolver::newMeasurementSeries()
    {
        iter = 0;

        for ( auto && measure : *convergenceMeasures )
            measure->newMeasurementSeries();

        fluid->couplingData.dataprev.setZero();
        solid->couplingData.dataprev.setZero();
        fluid->couplingData.data.setZero();
        solid->couplingData.data.setZero();
    }
}
