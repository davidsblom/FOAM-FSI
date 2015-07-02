
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ESDIRK.H"

namespace sdc
{
    ESDIRK::ESDIRK(
        std::shared_ptr<SDCSolver> solver,
        std::string method,
        std::shared_ptr<AdaptiveTimeStepper> adaptiveTimeStepper
        )
        :
        solver( solver ),
        adaptiveTimeStepper( adaptiveTimeStepper ),
        dt( solver->getTimeStep() ),
        nbStages( 0 ),
        A(),
        B(),
        C(),
        Bhat(),
        N( solver->getDOF() )
    {
        assert( method == "SDIRK2" || method == "SDIRK3" || method == "SDIRK4" || method == "ESDIRK3" || method == "ESDIRK4" || method == "ESDIRK5" );
        assert( solver );
        assert( adaptiveTimeStepper );
        assert( dt > 0 );

        // source: Ellsiepen. Habilitation thesis Philipp Birken
        if ( method == "SDIRK2" )
        {
            adaptiveTimeStepper->setOrderEmbeddedMethod( 1 );
            double alpha = 1.0 - 0.5 * std::sqrt( 2 );
            double alphahat = 2.0 - 5.0 / 4.0 * std::sqrt( 2 );
            A.resize( 2, 2 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = alpha;
            A( 1, 0 ) = 1.0 - alpha;
            A( 1, 1 ) = alpha;
            Bhat( 0 ) = 1 - alphahat;
            Bhat( 1 ) = alphahat;
        }

        // source: Cash. Habilitation thesis Philipp Birken
        if ( method == "SDIRK3" )
        {
            adaptiveTimeStepper->setOrderEmbeddedMethod( 2 );
            A.resize( 3, 3 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            double alpha = 1.2084966491760101;
            double beta = -0.6443631706844691;
            double gamma = 0.4358665215084580;
            double delta = 0.7179332607542295;
            double alphahat = 0.7726301276675511;
            double betahat = 0.2273698723324489;
            A( 0, 0 ) = gamma;
            A( 1, 0 ) = delta - gamma;
            A( 1, 1 ) = gamma;
            A( 2, 0 ) = alpha;
            A( 2, 1 ) = beta;
            A( 2, 2 ) = gamma;
            Bhat( 0 ) = alphahat;
            Bhat( 1 ) = betahat;
            Bhat( 2 ) = 0;
        }

        // source: Cash-5-3-4 http://runge.math.smu.edu/arkode_dev/doc/guide/build/html/Butcher.html
        if ( method == "SDIRK4" )
        {
            adaptiveTimeStepper->setOrderEmbeddedMethod( 3 );
            A.resize( 5, 5 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = 0.435866521508;
            A( 1, 0 ) = -1.13586652150;
            A( 1, 1 ) = 0.435866521508;
            A( 2, 0 ) = 1.08543330679;
            A( 2, 1 ) = -0.721299828287;
            A( 2, 2 ) = 0.435866521508;
            A( 3, 0 ) = 0.416349501547;
            A( 3, 1 ) = 0.190984004184;
            A( 3, 2 ) = -0.118643265417;
            A( 3, 3 ) = 0.435866521508;
            A( 4, 0 ) = 0.896869652944;
            A( 4, 1 ) = 0.0182725272734;
            A( 4, 2 ) = -0.0845900310706;
            A( 4, 3 ) = -0.266418670647;
            A( 4, 4 ) = 0.435866521508;
            Bhat( 0 ) = 0.776691932910;
            Bhat( 1 ) = 0.0297472791484;
            Bhat( 2 ) = -0.0267440239074;
            Bhat( 3 ) = 0.220304811849;
            Bhat( 4 ) = 0;
        }

        if ( method == "ESDIRK3" )
        {
            adaptiveTimeStepper->setOrderEmbeddedMethod( 2 );
            A.resize( 4, 4 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = 0.0;
            A( 1, 0 ) = 1767732205903.0 / 4055673282236.0;
            A( 1, 1 ) = 1767732205903.0 / 4055673282236.0;
            A( 2, 0 ) = 2746238789719.0 / 10658868560708.0;
            A( 2, 1 ) = -640167445237.0 / 6845629431997.0;
            A( 2, 2 ) = 1767732205903.0 / 4055673282236.0;
            A( 3, 0 ) = 1471266399579.0 / 7840856788654.0;
            A( 3, 1 ) = -4482444167858.0 / 7529755066697.0;
            A( 3, 2 ) = 11266239266428.0 / 11593286722821.0;
            A( 3, 3 ) = 1767732205903.0 / 4055673282236.0;
            Bhat( 0 ) = 2756255671327.e0 / 12835298489170.e0;
            Bhat( 1 ) = -10771552573575.e0 / 22201958757719.e0;
            Bhat( 2 ) = 9247589265047.e0 / 10645013368117.e0;
            Bhat( 3 ) = 2193209047091.e0 / 5459859503100.e0;
        }

        if ( method == "ESDIRK4" )
        {
            adaptiveTimeStepper->setOrderEmbeddedMethod( 3 );
            A.resize( 6, 6 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 1, 0 ) = 1.e0 / 4.e0;
            A( 1, 1 ) = 1.e0 / 4.e0;
            A( 2, 0 ) = 8611.e0 / 62500.e0;
            A( 2, 1 ) = -1743.e0 / 31250.e0;
            A( 2, 2 ) = 1.e0 / 4.e0;
            A( 3, 0 ) = 5012029.e0 / 34652500.e0;
            A( 3, 1 ) = -654441.e0 / 2922500.e0;
            A( 3, 2 ) = 174375.e0 / 388108.e0;
            A( 3, 3 ) = 1.e0 / 4.e0;
            A( 4, 0 ) = 15267082809.e0 / 155376265600.e0;
            A( 4, 1 ) = -71443401.e0 / 120774400.e0;
            A( 4, 2 ) = 730878875.e0 / 902184768.e0;
            A( 4, 3 ) = 2285395.e0 / 8070912.e0;
            A( 4, 4 ) = 1.e0 / 4.e0;
            A( 5, 0 ) = 82889.e0 / 524892.e0;
            A( 5, 1 ) = 0.e0;
            A( 5, 2 ) = 15625.e0 / 83664.e0;
            A( 5, 3 ) = 69875.e0 / 102672.e0;
            A( 5, 4 ) = -2260.e0 / 8211.e0;
            A( 5, 5 ) = 1.e0 / 4.e0;
            Bhat( 0 ) = 4586570599.e0 / 29645900160.e0;
            Bhat( 1 ) = 0.e0;
            Bhat( 2 ) = 178811875.e0 / 945068544.e0;
            Bhat( 3 ) = 814220225.e0 / 1159782912.e0;
            Bhat( 4 ) = -3700637.e0 / 11593932.e0;
            Bhat( 5 ) = 61727.e0 / 225920.e0;
        }

        if ( method == "ESDIRK5" )
        {
            adaptiveTimeStepper->setOrderEmbeddedMethod( 4 );
            A.resize( 8, 8 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = 0.e0;
            A( 1, 0 ) = 41.e0 / 200.e0;
            A( 1, 1 ) = 41.e0 / 200.e0;
            A( 2, 0 ) = 41.e0 / 400.e0;
            A( 2, 1 ) = -567603406766.e0 / 11931857230679.e0;
            A( 2, 2 ) = 41.e0 / 200.e0;
            A( 3, 0 ) = 683785636431.e0 / 9252920307686.e0;
            A( 3, 1 ) = 0.e0;
            A( 3, 2 ) = -110385047103.e0 / 1367015193373.e0;
            A( 3, 3 ) = 41.e0 / 200.e0;
            A( 4, 0 ) = 3016520224154.e0 / 10081342136671.e0;
            A( 4, 1 ) = 0.e0;
            A( 4, 2 ) = 30586259806659.e0 / 12414158314087.e0;
            A( 4, 3 ) = -22760509404356.e0 / 11113319521817.e0;
            A( 4, 4 ) = 41.e0 / 200.e0;
            A( 5, 0 ) = 218866479029.e0 / 1489978393911.e0;
            A( 5, 1 ) = 0.e0;
            A( 5, 2 ) = 638256894668.e0 / 5436446318841.e0;
            A( 5, 3 ) = -1179710474555.e0 / 5321154724896.e0;
            A( 5, 4 ) = -60928119172.e0 / 8023461067671.e0;
            A( 5, 5 ) = 41.e0 / 200.e0;
            A( 6, 0 ) = 1020004230633.e0 / 5715676835656.e0;
            A( 6, 1 ) = 0.e0;
            A( 6, 2 ) = 25762820946817.e0 / 25263940353407.e0;
            A( 6, 3 ) = -2161375909145.e0 / 9755907335909.e0;
            A( 6, 4 ) = -211217309593.e0 / 5846859502534.e0;
            A( 6, 5 ) = -4269925059573.e0 / 7827059040749.e0;
            A( 6, 6 ) = 41.e0 / 200.e0;
            A( 7, 0 ) = -872700587467.e0 / 9133579230613.e0;
            A( 7, 1 ) = 0.e0;
            A( 7, 2 ) = 0.e0;
            A( 7, 3 ) = 22348218063261.e0 / 9555858737531.e0;
            A( 7, 4 ) = -1143369518992.e0 / 8141816002931.e0;
            A( 7, 5 ) = -39379526789629.e0 / 19018526304540.e0;
            A( 7, 6 ) = 32727382324388.e0 / 42900044865799.e0;
            A( 7, 7 ) = 41.e0 / 200.e0;
            Bhat( 0 ) = -975461918565.e0 / 9796059967033.e0;
            Bhat( 1 ) = 0.e0;
            Bhat( 2 ) = 0.e0;
            Bhat( 3 ) = 78070527104295.e0 / 32432590147079.e0;
            Bhat( 4 ) = -548382580838.e0 / 3424219808633.e0;
            Bhat( 5 ) = -33438840321285.e0 / 15594753105479.e0;
            Bhat( 6 ) = 3629800801594.e0 / 4656183773603.e0;
            Bhat( 7 ) = 4035322873751.e0 / 18575991585200.e0;
        }

        nbStages = A.cols();

        C.resize( nbStages );

        for ( int i = 0; i < nbStages; i++ )
            C( i ) = A.row( i ).sum();

        B.resize( nbStages );

        for ( int i = 0; i < nbStages; i++ )
            B( i ) = A( nbStages - 1, i );

        solver->setNumberOfStages( nbStages );

        adaptiveTimeStepper->setEndTime( solver->getEndTime() );
    }

    ESDIRK::~ESDIRK()
    {}

    bool ESDIRK::isStageImplicit( double Akk )
    {
        return Akk > 0;
    }

    void ESDIRK::run()
    {
        double t = 0;

        while ( std::abs( t - solver->getEndTime() ) > 1.0e-13 && t < solver->getEndTime() )
        {
            double computedTimeStep = dt;

            solveTimeStep( t );

            if ( adaptiveTimeStepper->isAccepted() )
                t += computedTimeStep;
        }
    }

    void ESDIRK::solveTimeStep( const double t0 )
    {
        if ( adaptiveTimeStepper->isPreviousStepAccepted() )
            solver->nextTimeStep();

        Eigen::MatrixXd solStages( nbStages, N ), F( nbStages, N );
        F.setZero();

        Eigen::VectorXd sol( N ), f( N ), qold( N ), result( N ), rhs( N );
        solver->getSolution( sol );
        solStages.row( 0 ) = sol;

        qold = sol;

        double t = t0;

        solver->evaluateFunction( 0, sol, t, f );
        F.row( 0 ) = f;

        // Keep the solution of the first stage and function evaluate
        // in memory in case the time step is rejected
        Eigen::VectorXd solOld = sol;
        Eigen::VectorXd fOld = f;

        // Loop over the stages

        solver->initTimeStep();

        for ( int j = 0; j < nbStages; j++ )
        {
            if ( !isStageImplicit( A( j, j ) ) )
                continue;

            t = t0 + C( j ) * dt;

            Info << "\nTime = " << t << ", ESDIRK stage = " << j + 1 << "/" << nbStages << nl << endl;

            rhs.setZero();

            // Calculate sum of the stage residuals
            for ( int iStage = 0; iStage < j; iStage++ )
                rhs += A( j, iStage ) * F.row( iStage ).transpose();

            rhs.array() *= dt;

            solver->implicitSolve( false, j, t, A( j, j ) * dt, qold, rhs, f, result );

            solStages.row( j ) = result;
            F.row( j ) = f;
        }

        if ( adaptiveTimeStepper->isEnabled() )
        {
            double newTimeStep = 0;
            Eigen::VectorXd errorEstimate( N );
            errorEstimate.setZero();

            for ( int iStage = 0; iStage < nbStages; iStage++ )
                errorEstimate += ( B( iStage ) - Bhat( iStage ) ) * F.row( iStage );

            errorEstimate *= dt;

            bool accepted = adaptiveTimeStepper->determineNewTimeStep( errorEstimate, result, dt, newTimeStep );

            dt = newTimeStep;

            if ( not accepted )
                solver->setSolution( solOld, fOld );
        }

        if ( adaptiveTimeStepper->isAccepted() )
            solver->finalizeTimeStep();
    }
}
