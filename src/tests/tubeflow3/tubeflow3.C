
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <fstream>
#include <iostream>
#include <deque>
#include <string>

#include "ASMILS.H"
#include "AggressiveSpaceMapping.H"
#include "BroydenPostProcessing.H"
#include "MLIQNILSSolver.H"
#include "ManifoldMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MultiLevelSpaceMappingSolver.H"
#include "OutputSpaceMapping.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "AbsoluteConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "TubeFlowLinearizedSolidSolver.H"
#include "TubeFlowLinearizedFluidSolver.H"
#include "AndersonPostProcessing.H"

using namespace std;
using namespace fsi;
using namespace tubeflow;

int main()
{
    // List of solvers
    deque<std::string> fsiSolvers;

    fsiSolvers.push_back( "Anderson" );
    fsiSolvers.push_back( "MM" );

    std::string fsiSolver = fsiSolvers.at( 1 );

    // Physical settings
    scalar r0 = 0.2;
    scalar a0 = M_PI * r0 * r0;
    scalar u0 = 0.1;
    scalar p0 = 0;
    scalar dt = 0.1;
    scalar L = 1;
    scalar T = 10;
    scalar rho_f = 1.225;
    scalar rho_s = 1.225;
    scalar nu = 0.5;
    scalar E = 490;
    scalar G = 490;
    scalar h = 1.0e-3;
    scalar cmk = std::sqrt( E * h / (2 * rho_f * r0) );

    // Computational settings
    scalar tol = 1.0e-5;
    int maxIter = 100;
    scalar initialRelaxation = 1.0e-3;
    scalar singularityLimit = 1.0e-14;
    int extrapolation = 2;
    int minIter = 1;
    int reuseInformationStartingFromTimeIndex = 0;
    bool scaling = false;
    scalar beta = 1;
    bool updateJacobian = false;

    // Grid settings
    int N = 1000;
    int couplingGridSize = N;

    int nbReuse = 4;
    bool parallel = false;
    int nbLevels = 2;

    std::string label = fsiSolver + "_" + to_string( nbReuse );
    label += "_" + to_string( nbLevels );

    ifstream ifile( label + ".log" );

    // Create shared pointers to solvers
    shared_ptr<BaseMultiLevelSolver> fluid;
    shared_ptr<BaseMultiLevelSolver> solid;
    shared_ptr<BaseMultiLevelSolver> fineModelFluid;
    shared_ptr<BaseMultiLevelSolver> fineModelSolid;
    shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
    shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
    shared_ptr<MultiLevelFsiSolver> fsi;
    shared_ptr< list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
    shared_ptr<PostProcessing> postProcessing;
    shared_ptr< deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;
    shared_ptr< deque<shared_ptr<SpaceMappingSolver> > > solvers;
    shared_ptr<Solver> solver;
    shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

    models = shared_ptr< deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > ( new deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );
    solvers = shared_ptr< deque<shared_ptr<SpaceMappingSolver> > > ( new deque<shared_ptr<SpaceMappingSolver> >() );

    if ( fsiSolver == "MM" )
    {
        tol = 1.0e-5;
        N = couplingGridSize;

        fineModelFluid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
        fineModelSolid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E, r0, T ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, nbLevels - 1 ) );
        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, 1, nbLevels - 1 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        int maxUsedIterations = fineModelSolid->data.rows() * fineModelSolid->data.cols();

        postProcessing = shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        solid.reset();
        multiLevelFluidSolver.reset();
        multiLevelSolidSolver.reset();
        convergenceMeasures.reset();
        multiLevelFsiSolver.reset();
        postProcessing.reset();

        for ( int level = 0; level < nbLevels - 1; level++ )
        {
            if ( level == 0 )
            {
                tol = 1.0e-1;
                N = 100;
            }

            if ( level == 1 )
            {
                tol = 1.0e-5;
                N = couplingGridSize;
            }

            fluid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowLinearizedFluidSolver( N, p0, r0, u0, rho_f, E, h, T, dt, L ) );
            solid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E, r0, T ) );

            multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, 0, level ) );
            multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, 1, level ) );

            // Convergence measures
            convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new ResidualRelativeConvergenceMeasure( 0, true, tol ) ) );
            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new AbsoluteConvergenceMeasure( 0, true, 1.0e-17 ) ) );

            multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

            int maxUsedIterations = solid->data.rows() * solid->data.cols();

            postProcessing = shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

            shared_ptr<ImplicitMultiLevelFsiSolver> implicitMultiLevelFsiSolver( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

            models->push_back( implicitMultiLevelFsiSolver );

            solid.reset();
            multiLevelFluidSolver.reset();
            multiLevelSolidSolver.reset();
            convergenceMeasures.reset();
            multiLevelFsiSolver.reset();
            postProcessing.reset();
            implicitMultiLevelFsiSolver.reset();
        }

        models->push_back( fineModel );

        assert( static_cast<int>( models->size() ) == nbLevels );

        if ( fsiSolver == "MM" )
        {
            for ( int level = 0; level < nbLevels - 1; level++ )
            {
                shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
                shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel;

                fineModel = models->at( level + 1 );

                if ( level == 0 )
                    coarseModel = models->at( level );

                if ( level > 0 )
                    coarseModel = solvers->at( level - 1 );

                shared_ptr<SpaceMapping> spaceMapping;

                spaceMapping = shared_ptr<SpaceMapping> ( new ManifoldMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian ) );

                shared_ptr<SpaceMappingSolver > spaceMappingSolver( new SpaceMappingSolver( fineModel, coarseModel, spaceMapping ) );

                solvers->push_back( spaceMappingSolver );

                spaceMapping.reset();
                spaceMappingSolver.reset();
            }

            assert( static_cast<int>( solvers->size() ) == nbLevels - 1 );
        }
    }

    if ( fsiSolver == "Anderson" )
    {
        assert( !fsi );

        tol = 1.0e-5;
        N = couplingGridSize;

        fluid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
        // solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
        solid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E, r0, T ) );

        shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, 0, 0 ) );
        shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, 1, 0 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        int maxUsedIterations = solid->data.rows() * solid->data.cols();

        fsi = shared_ptr<MultiLevelFsiSolver>( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<PostProcessing> ( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        solver = shared_ptr<Solver>( new ImplicitMultiLevelFsiSolver( fsi, postProcessing ) );
    }

    if ( fsiSolver == "MM" )
        solver = shared_ptr<Solver> ( new MultiLevelSpaceMappingSolver( solvers, models, true ) );

    try
    {
        solver->run();
    }
    catch ( ... )
    {

    }
    ofstream logFile( label + ".log" );
    assert( logFile.is_open() );

    logFile << "label = " << label << endl;
    logFile << "nbLevels = " << nbLevels << endl;
    logFile << "solver = " << fsiSolver << endl;
    logFile << "nbReuse = " << nbReuse << endl;
    logFile << "parallel = " << parallel << endl;
    logFile << "elapsed time = " << solver->timeElapsed() << endl;
    logFile << "nbTimeSteps = " << fluid->timeIndex << endl;
    logFile << "updateJacobian = " << updateJacobian << endl;

    if ( fsiSolver == "Anderson" )
        logFile << "nbIter = " << fsi->nbIter << endl;

    if ( fsiSolver == "MM" )
    {
        int level = 0;

        for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = models->begin(); it != models->end(); ++it )
        {
            shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;
            logFile << "level " << level << " nbIter = " << model->fsi->nbIter << endl;
            logFile << "level " << level << " N = " << model->fsi->fluid->data.rows() << endl;

            level++;
        }
    }

    logFile.close();

    return 0;
}
