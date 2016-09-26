
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "BroydenPostProcessing.H"
#include "ManifoldMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "MultiLevelFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "SpaceMappingSolver.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "AndersonPostProcessing.H"
#include "gtest/gtest.h"
#include "LinearFunction.H"
#include <iomanip>

using namespace tubeflow;
using namespace rbf;
using std::shared_ptr;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class ManifoldMappingSolverTest2 : public::testing::Test
{
protected:

    virtual void SetUp()
    {
        // Physical settings
        scalar r0 = 0.2;
        std::cout.precision( 15 );
        std::cout << "M_PI = " << M_PI << std::endl;
        scalar a0 = M_PI * r0 * r0;
        std::cout << "a0 = " << a0 << std::endl;
        scalar u0 = 0.1;
        scalar p0 = 0;
        scalar dt = 0.1;
        int N = 10;
        scalar L = 1;
        scalar T = 10;
        scalar rho = 1.225;
        scalar E = 490;
        scalar h = 1.0e-3;
        scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
        scalar c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
        scalar kappa = c0 / u0;
        scalar tau = u0 * dt / L;

        // Computational settings
        scalar tol = 1.0e-4;
        int maxIter = 50;
        scalar initialRelaxation = 1.0e-3;
        scalar singularityLimit = 1.0e-17;
        int reuseInformationStartingFromTimeIndex = 0;
        scalar beta = 1;
        bool updateJacobian = false;

        // Parametrized settings
        bool parallel = true;
        int nbReuse = 0;
        int extrapolation = 0;
        int couplingGridSize = 100;
        bool scaling = false;

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );

        // Initialize variables
        shared_ptr<TubeFlowSolidSolver> solid;
        shared_ptr<TubeFlowFluidSolver> fluid;
        shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        shared_ptr<PostProcessing> postProcessing;
        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

        // Fine model

        int maxUsedIterations = couplingGridSize;

        fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, 1 ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, solid, 1, 1 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Coarse model

        maxUsedIterations = N;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, 0, 0 ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->solid, 1, 0 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 0.1 * tol ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 0.1 * tol ) ) );
        // convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure(1,false, 20)) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        shared_ptr<ManifoldMapping> manifoldMapping( new ManifoldMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian, true ) );

        // Create manifold mapping solver
        solver = new SpaceMappingSolver( fineModel, coarseModel, manifoldMapping );

        assert( fluid );
        assert( solid );
        assert( multiLevelFluidSolver );
        assert( multiLevelSolidSolver );
        assert( convergenceMeasures );
        assert( multiLevelFsiSolver );
        assert( postProcessing );
    }

    virtual void TearDown()
    {
        delete solver;
    }

    SpaceMappingSolver * solver;
    shared_ptr<TubeFlowFluidSolver> fineModelFluid;
};

TEST_F( ManifoldMappingSolverTest2, object )
{
    ASSERT_TRUE( true );
}

TEST_F( ManifoldMappingSolverTest2, solveTimeStep )
{
    for ( int i = 0; i < 100; i++)
    {
        solver->solveTimeStep();
        Info << "Fine model iterations = " << solver->fineModel->fsi->nbIter << endl;
        Info << "Coarse model iterations = " << solver->coarseModel->fsi->nbIter - solver->fineModel->fsi->nbIter << endl;
    }
}
