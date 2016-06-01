
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

    // Physical settings
    const scalar r0 = 0.2;
    const scalar a0 = M_PI * r0 * r0;
    const scalar u0 = 0.1;
    const scalar p0 = 0;
    const scalar dt = 0.1;
    const scalar L = 1;
    const scalar T = 10;
    const scalar rho_f = 1.225;
    const scalar rho_s = 1.225;
    const scalar E = 490;
    const scalar h = 1.0e-3;
    const scalar cmk = std::sqrt( E * h / (2 * rho_f * r0) );
    const scalar G = 490;
    const scalar nu = 0.5;

    // Computational settings
    const scalar initialRelaxation = 1.0e-3;
    const scalar singularityLimit = 1.0e-14;
    const int extrapolation = 0;
    const int minIter = 1;
    const int reuseInformationStartingFromTimeIndex = 0;
    const bool scaling = false;
    const scalar beta = 0.5;
    const scalar tol = 1.0e-5;
    const int maxIter = 20;
    const bool parallel = false;
    const bool initialSolutionCoarseModel = true;

    std::deque<std::deque<int> > levelsList = { {100}, {20, 100}, {50, 100}, {80, 100}, {100, 100}, {50, 80, 100}, {20, 50, 100}, {20, 80, 100} };

    for ( std::string fluidSolverSetting : {"linearized", "non-linear"} )
    {
        for ( auto updateJacobian : {true, false} )
        {
            for ( auto nbReuse : {0, 1, 4} )
            {
                for ( auto fsiSolver : fsiSolvers )
                {
                    for ( auto levels : levelsList )
                    {
                        const int couplingGridSize = levels.back();
                        int N = couplingGridSize;
                        const int nbLevels = levels.size();

                        if ( nbLevels > 1 && fsiSolver == "Anderson" )
                            continue;

                        if ( fluidSolverSetting != "non-linear" && fsiSolver == "Anderson" )
                            continue;

                        if ( nbLevels == 1 && fsiSolver == "MM" )
                            continue;

                        if ( updateJacobian && nbReuse > 0 )
                            continue;

                        std::string label = fsiSolver;

                        if ( updateJacobian )
                            label += "_complete_jacobian";
                        else
                            label += "_reuse_" + to_string( nbReuse );

                        if ( fsiSolver == "MM" )
                        {
                            label += "_" + to_string( nbLevels );

                            for ( auto level : levels )
                                label += "_" + to_string( level );

                            label += "_fluid_" + fluidSolverSetting;
                        }

                        ifstream ifile( label + ".log" );

                        if ( ifile )
                            continue; // log file exists;

                        // Create shared pointers to solvers

                        shared_ptr<BaseMultiLevelSolver> fineModelFluid;
                        shared_ptr<BaseMultiLevelSolver> fineModelSolid;

                        shared_ptr<MultiLevelFsiSolver> fsi;
                        shared_ptr< list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
                        shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
                        shared_ptr< deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;
                        shared_ptr< deque<shared_ptr<SpaceMappingSolver> > > solvers;
                        shared_ptr<Solver> solver;
                        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

                        models = shared_ptr< deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > ( new deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );
                        solvers = shared_ptr< deque<shared_ptr<SpaceMappingSolver> > > ( new deque<shared_ptr<SpaceMappingSolver> >() );

                        fineModelFluid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
                        fineModelSolid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E, r0, T ) );

                        // fineModelSolid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho_s, L, N ) );

                        if ( fsiSolver == "MM" )
                        {
                            shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
                            shared_ptr<MultiLevelSolver> multiLevelSolidSolver;

                            multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, nbLevels - 1 ) );
                            multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, 1, nbLevels - 1 ) );

                            // Convergence measures
                            convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

                            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
                            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

                            multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

                            int maxUsedIterations = fineModelSolid->data.rows() * fineModelSolid->data.cols();

                            shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

                            fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

                            multiLevelFluidSolver.reset();
                            multiLevelSolidSolver.reset();
                            convergenceMeasures.reset();
                            multiLevelFsiSolver.reset();
                            postProcessing.reset();

                            for ( int level = 0; level < nbLevels - 1; level++ )
                            {
                                N = levels[level];

                                shared_ptr<BaseMultiLevelSolver> fluid;

                                if ( fluidSolverSetting == "non-linear" )
                                    fluid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );

                                if ( fluidSolverSetting == "linearized" )
                                    fluid = shared_ptr<BaseMultiLevelSolver> ( new TubeFlowLinearizedFluidSolver( N, p0, r0, u0, rho_f, E, h, T, dt, L ) );

                                shared_ptr<BaseMultiLevelSolver> solid( new TubeFlowLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E, r0, T ) );

                                shared_ptr<MultiLevelSolver> multiLevelFluidSolver( new MultiLevelSolver( fluid, fineModelFluid, 0, level ) );
                                shared_ptr<MultiLevelSolver> multiLevelSolidSolver( new MultiLevelSolver( solid, fineModelSolid, 1, level ) );

                                // Convergence measures
                                convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

                                convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
                                convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new ResidualRelativeConvergenceMeasure( 0, true, 1.0e-2 ) ) );
                                convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new AbsoluteConvergenceMeasure( 0, true, 1.0e-17 ) ) );

                                multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

                                int maxUsedIterations = solid->data.rows() * solid->data.cols();

                                shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

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

                                    spaceMapping = shared_ptr<SpaceMapping> ( new ManifoldMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian, initialSolutionCoarseModel ) );

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

                            shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, 0 ) );
                            shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( fineModelSolid, fineModelFluid, 1, 0 ) );

                            // Convergence measures
                            convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

                            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
                            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

                            int maxUsedIterations = fineModelSolid->data.rows() * fineModelSolid->data.cols();

                            fsi = shared_ptr<MultiLevelFsiSolver>( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

                            shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

                            solver = shared_ptr<Solver>( new ImplicitMultiLevelFsiSolver( fsi, postProcessing ) );
                        }

                        if ( fsiSolver == "MM" )
                            solver = shared_ptr<Solver> ( new MultiLevelSpaceMappingSolver( solvers, models, true ) );

                        std::cout << "label = " << label << std::endl;

                        solver->run();

                        ofstream logFile( label + ".log" );
                        assert( logFile.is_open() );

                        logFile << "label = " << label << endl;
                        logFile << "solver = " << fsiSolver << endl;
                        logFile << "nbReuse = " << nbReuse << endl;
                        logFile << "parallel = " << parallel << endl;
                        logFile << "elapsed time = " << solver->timeElapsed() << endl;
                        logFile << "nbTimeSteps = " << fineModelFluid->timeIndex << endl;
                        logFile << "updateJacobian = " << updateJacobian << endl;

                        if ( fsiSolver == "Anderson" )
                        {
                            logFile << "nbIter = " << fsi->nbIter << endl;
                            logFile << "avgIter = " << scalar( fsi->nbIter ) / scalar( fineModelFluid->timeIndex ) << endl;
                        }

                        if ( fsiSolver == "MM" )
                        {
                            logFile << "nbLevels = " << nbLevels << endl;

                            int level = 0;

                            for ( auto && model : *models )
                            {
                                logFile << "level " << level << " nbIter = " << model->fsi->nbIter << endl;
                                logFile << "level " << level << " avgIter = " << scalar( model->fsi->nbIter ) / scalar( fineModelFluid->timeIndex ) << endl;
                                logFile << "level " << level << " N = " << model->fsi->fluid->data.rows() << endl;

                                level++;
                            }
                        }

                        logFile.close();
                    }
                }
            }
        }
    }

    return 0;
}
