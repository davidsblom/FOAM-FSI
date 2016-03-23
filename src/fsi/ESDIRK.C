
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
        N( solver->getDOF() ),
        stageIndex( 0 ),
        F(),
        solStages(),
        qold(),
        timeIndex( 0 )
    {
        assert( solver );
        assert( adaptiveTimeStepper );
        assert( dt > 0 );

        initializeButcherTableau( method );

        bool firstStageImplicit = isStageImplicit( A( 0, 0 ) );

        if ( firstStageImplicit )
            solver->setNumberOfImplicitStages( nbStages );
        else
            solver->setNumberOfImplicitStages( nbStages - 1 );

        adaptiveTimeStepper->setEndTime( solver->getEndTime() );
    }

    ESDIRK::ESDIRK( std::string method )
        :
        solver( nullptr ),
        adaptiveTimeStepper( nullptr ),
        dt( -1 ),
        nbStages( 0 ),
        A(),
        B(),
        C(),
        Bhat(),
        N( 0 ),
        stageIndex( 0 ),
        F(),
        solStages(),
        qold(),
        timeIndex( 0 )
    {
        initializeButcherTableau( method );
    }

    ESDIRK::~ESDIRK()
    {}

    bool ESDIRK::isStageImplicit( scalar Akk )
    {
        return Akk > 0;
    }

    void ESDIRK::run()
    {
        scalar t = solver->getStartTime();

        while ( std::abs( t - solver->getEndTime() ) > 1.0e-13 && t < solver->getEndTime() )
        {
            scalar computedTimeStep = dt;

            solveTimeStep( t );

            if ( adaptiveTimeStepper->isAccepted() )
                t += computedTimeStep;
        }
    }

    void ESDIRK::solveTimeStep( const scalar t0 )
    {
        if ( adaptiveTimeStepper->isPreviousStepAccepted() )
            solver->nextTimeStep();

        fsi::matrix solStages( nbStages, N ), F( nbStages, N );
        F.setZero();

        fsi::vector sol( N ), f( N ), qold( N ), result( N ), rhs( N );
        solver->getSolution( sol, f );
        solStages.row( 0 ) = sol;

        qold = sol;

        scalar t = t0;

        solver->evaluateFunction( 0, sol, t, f );
        F.row( 0 ) = f;

        // Keep the solution of the first stage and function evaluate
        // in memory in case the time step is rejected
        fsi::vector solOld = sol;
        fsi::vector fOld = f;

        // Loop over the stages

        solver->initTimeStep();

        int iImplicitStage = 0;

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

            solver->implicitSolve( false, iImplicitStage, 0, t, A( j, j ) * dt, qold, rhs, f, result );

            assert( (1.0 / (A( j, j ) * dt) * (result - qold - rhs) - f).norm() < 1.0e-12 );

            solStages.row( j ) = result;
            F.row( j ) = f;

            iImplicitStage++;
        }

        if ( adaptiveTimeStepper->isEnabled() )
        {
            scalar newTimeStep = 0;
            fsi::vector errorEstimate( N );
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

    void ESDIRK::getSourceTerm(
        const bool /*corrector*/,
        const int k,
        const int /*sweep*/,
        const scalar deltaT,
        fsi::vector & rhs,
        fsi::vector & qold
        )
    {
        assert( k <= nbStages - 1 );

        qold = this->qold;

        // Compute the time step from the stage deltaT
        if ( dt < 0 )
        {
            // first time step, first prediction step
            dt = deltaT / A( k, k );
        }

        assert( dt > 0 );

        rhs.setZero();

        // Calculate sum of the stage residuals
        for ( int iStage = 0; iStage < k; iStage++ )
            rhs += A( k, iStage ) * F.row( iStage ).transpose();

        rhs.array() *= dt;

        this->stageIndex = k;
    }

    void ESDIRK::setFunction(
        const int k,
        const fsi::vector & f,
        const fsi::vector & result
        )
    {
        assert( f.rows() == result.rows() );
        assert( k <= nbStages - 1 );

        if ( F.cols() == 0 )
            F.resize( nbStages, f.rows() );

        if ( solStages.cols() == 0 )
            solStages.resize( nbStages, f.rows() );

        F.row( k ) = f;
        solStages.row( k ) = result;
    }

    void ESDIRK::setOldSolution(
        int timeIndex,
        const fsi::vector & result
        )
    {
        assert( timeIndex >= this->timeIndex );

        if ( qold.rows() == result.rows() )
            qold = solStages.bottomRows( 1 ).transpose();
        else
        if ( timeIndex > this->timeIndex )
            qold = result;

        this->timeIndex = timeIndex;
    }

    int ESDIRK::getNbImplicitStages()
    {
        assert( nbStages > 0 );

        bool firstStageImplicit = isStageImplicit( A( 0, 0 ) );

        if ( firstStageImplicit )
            return nbStages;
        else
            return nbStages - 1;
    }

    void ESDIRK::initializeButcherTableau( std::string method )
    {
        assert( method == "SDIRK2" || method == "SDIRK3" || method == "SDIRK4" || method == "ESDIRK3" || method == "ESDIRK4" || method == "ESDIRK5" || method == "ESDIRK53PR" || method == "ESDIRK63PR" || method == "ESDIRK74PR" );

        // source: Ellsiepen. Habilitation thesis Philipp Birken
        if ( method == "SDIRK2" )
        {
            if ( adaptiveTimeStepper )
                adaptiveTimeStepper->setOrderEmbeddedMethod( 1 );

            scalar alpha = 1.0 - 0.5 * std::sqrt( 2 );
            scalar alphahat = 2.0 - 5.0 / 4.0 * std::sqrt( 2 );
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
            if ( adaptiveTimeStepper )
                adaptiveTimeStepper->setOrderEmbeddedMethod( 2 );

            A.resize( 3, 3 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            scalar alpha = 1.2084966491760101;
            scalar beta = -0.6443631706844691;
            scalar gamma = 0.4358665215084580;
            scalar delta = 0.7179332607542295;
            scalar alphahat = 0.7726301276675511;
            scalar betahat = 0.2273698723324489;
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
            if ( adaptiveTimeStepper )
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
            if ( adaptiveTimeStepper )
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
            if ( adaptiveTimeStepper )
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
            if ( adaptiveTimeStepper )
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

        // source: http://www.sciencedirect.com/science/article/pii/S0168927415000501
        if ( method == "ESDIRK53PR" )
        {
            if ( adaptiveTimeStepper )
                adaptiveTimeStepper->setOrderEmbeddedMethod( 2 );

            A.resize( 5, 5 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = 0.0;
            A( 1, 0 ) = 2.777777777777778e-01;
            A( 1, 1 ) = 2.777777777777778e-01;
            A( 2, 0 ) = 3.456552483519272e-01;
            A( 2, 1 ) = 1.681740315717733e-01;
            A( 2, 2 ) = 2.777777777777778e-01;
            A( 3, 0 ) = 3.965643047257401e-01;
            A( 3, 1 ) = 1.001154404932533e-01;
            A( 3, 2 ) = 1.255424770032288e-01;
            A( 3, 3 ) = 2.777777777777778e-01;
            A( 4, 0 ) = 2.481479828780141e-01;
            A( 4, 1 ) = 2.139473588935955e-01;
            A( 4, 2 ) = 1.206274239267400e+00;
            A( 4, 3 ) = -9.461473588167871e-01;
            A( 4, 4 ) = 2.777777777777778e-01;
            Bhat( 0 ) = 4.445537532713554e-01;
            Bhat( 1 ) = -1.065203443758999e-01;
            Bhat( 2 ) = 2.533129069755295e-01;
            Bhat( 3 ) = 5.000000000000000e-01;
            Bhat( 4 ) = -9.134631587098500e-02;
        }

        // source: http://www.sciencedirect.com/science/article/pii/S0168927415000501
        if ( method == "ESDIRK63PR" )
        {
            if ( adaptiveTimeStepper )
                adaptiveTimeStepper->setOrderEmbeddedMethod( 2 );

            A.resize( 6, 6 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = 0.0;
            A( 1, 0 ) = 4.166666666666667e-01;
            A( 1, 1 ) = 4.166666666666667e-01;
            A( 2, 0 ) = 3.640473915723038e-01;
            A( 2, 1 ) = -4.189886135331312e-02;
            A( 2, 2 ) = 4.166666666666667e-01;
            A( 3, 0 ) = -2.894969214392781e+00;
            A( 3, 1 ) = -2.256341718064659e+01;
            A( 3, 2 ) = 2.534171972837271e+01;
            A( 3, 3 ) = 4.166666666666667e-01;
            A( 4, 0 ) = 2.309551022782098e-01;
            A( 4, 1 ) = -1.849667242832423e+00;
            A( 4, 2 ) = 2.197073089164931e+00;
            A( 4, 3 ) = 4.972384722615363e-03;
            A( 4, 4 ) = 4.166666666666667e-01;
            A( 5, 0 ) = 3.054968378466108e-01;
            A( 5, 1 ) = 4.057983152922798e+00;
            A( 5, 2 ) = -2.202162095667910e+00;
            A( 5, 3 ) = 1.333484429273537e-01;
            A( 5, 4 ) = -1.711333004695519e+00;
            A( 5, 5 ) = 4.166666666666667e-01;
            Bhat( 0 ) = 2.309551022782098e-01;
            Bhat( 1 ) = -1.849667242832423e+00;
            Bhat( 2 ) = 2.197073089164931e+00;
            Bhat( 3 ) = 4.972384722615363e-03;
            Bhat( 4 ) = 4.166666666666667e-01;
            Bhat( 5 ) = 0.0;
        }

        // source: http://www.sciencedirect.com/science/article/pii/S0168927415000501
        if ( method == "ESDIRK74PR" )
        {
            if ( adaptiveTimeStepper )
                adaptiveTimeStepper->setOrderEmbeddedMethod( 3 );

            A.resize( 7, 7 );
            A.setZero();
            Bhat.resize( A.cols() );
            Bhat.setZero();
            A( 0, 0 ) = 0.0;
            A( 1, 0 ) = 1.666666666666667e-01;
            A( 1, 1 ) = 1.666666666666667e-01;
            A( 2, 0 ) = 4.166666666666666e-02;
            A( 2, 1 ) = -4.166666666666666e-02;
            A( 2, 2 ) = 1.666666666666667e-01;
            A( 3, 0 ) = -1.500000000000000e+00;
            A( 3, 1 ) = -1.333333333333333e+00;
            A( 3, 2 ) = 3.333333333333333e+00;
            A( 3, 3 ) = 1.666666666666667e-01;
            A( 4, 0 ) = -1.580729166666667e+00;
            A( 4, 1 ) = -1.349609375000000e+00;
            A( 4, 2 ) = 3.472656250000000e+00;
            A( 4, 3 ) = 4.101562500000000e-02;
            A( 4, 4 ) = 1.666666666666667e-01;
            A( 5, 0 ) = -2.005366150605651e+00;
            A( 5, 1 ) = -1.768688648609954e+00;
            A( 5, 2 ) = 4.341269295345690e+00;
            A( 5, 3 ) = 2.326169434610579e-02;
            A( 5, 4 ) = 1.000000000000000e-01;
            A( 5, 5 ) = 1.666666666666667e-01;
            A( 6, 0 ) = 1.684854267805816e-01;
            A( 6, 1 ) = 7.501080898831836e-01;
            A( 6, 2 ) = -2.255843889686931e-01;
            A( 6, 3 ) = -9.134421504267402e-01;
            A( 6, 4 ) = 1.618140253772232e+00;
            A( 6, 5 ) = -5.643738977072310e-01;
            A( 6, 6 ) = 1.666666666666667e-01;
            Bhat( 0 ) = -3.930182461751728e-01;
            Bhat( 1 ) = 1.000000000000000e-01;
            Bhat( 2 ) = 9.916346405575472e-01;
            Bhat( 3 ) = 0.0;
            Bhat( 4 ) = -2.511232158528943e-01;
            Bhat( 5 ) = 4.393912810497486e-01;
            Bhat( 6 ) = 1.131155404207712e-01;
        }

        nbStages = A.cols();

        C.resize( nbStages );

        for ( int i = 0; i < nbStages; i++ )
            C( i ) = A.row( i ).sum();

        B.resize( nbStages );

        for ( int i = 0; i < nbStages; i++ )
            B( i ) = A( nbStages - 1, i );
    }

    void ESDIRK::outputResidual( std::string )
    {}

    bool ESDIRK::isConverged()
    {
        return true;
    }
}
