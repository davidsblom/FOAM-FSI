
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AbsoluteConvergenceMeasure.H"
#include "AndersonPostProcessing.H"
#include "ManifoldMapping.H"
#include "OutputSpaceMapping.H"
#include "AggressiveSpaceMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MultiLevelFsiSolver.H"
#include "MultiLevelSpaceMappingSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "SpaceMappingSolver.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;
using fsi::AndersonPostProcessing;
using fsi::MultiLevelSolver;
using fsi::MultiLevelFsiSolver;
using fsi::ImplicitMultiLevelFsiSolver;
using fsi::ConvergenceMeasure;
using fsi::PostProcessing;
using fsi::ManifoldMapping;
using fsi::SpaceMapping;
using fsi::SpaceMappingSolver;
using fsi::OutputSpaceMapping;
using fsi::AggressiveSpaceMapping;
using fsi::MinIterationConvergenceMeasure;
using fsi::RelativeConvergenceMeasure;

class MultiLevelSpaceMappingSolverParametrizedTest : public TestWithParam< std::tr1::tuple<int, int, int, int> >
{
protected:

    virtual void SetUp()
    {
        // Physical settings
        scalar r0 = 0.2;
        scalar a0 = M_PI * r0 * r0;
        scalar u0 = 0.1;
        scalar p0 = 0;
        scalar dt = 0.1;
        int N = 15;
        scalar L = 1;
        scalar T = 10;
        scalar dx = L / N;
        scalar rho = 1.225;
        scalar E = 490;
        scalar h = 1.0e-3;
        scalar cmk = std::sqrt( E * h / (2 * rho * r0) );
        scalar c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
        scalar kappa = c0 / u0;
        scalar tau = u0 * dt / L;

        // Computational settings
        scalar tol = 1.0e-4;
        scalar tolLiveSelection = 1.0e-4;
        int maxIter = 50;
        scalar initialRelaxation = 1.0e-3;
        scalar singularityLimit = 1.0e-14;
        bool scaling = false;
        bool updateJacobian = false;
        scalar beta = 0.5;
        int coarseningMinPoints = 5;
        int coarseningMaxPoints = 2000;
        bool parallel = false;
        int order = 0;
        int minIter = 5;
        int nbReuse = 2;
        bool coarsening = false;
        bool liveSelection = false;

        // Parametrized settings
        int extrapolation = std::tr1::get<0>( GetParam() );
        int couplingGridSize = std::tr1::get<1>( GetParam() );
        int spaceMappingAlgorithm = std::tr1::get<2>( GetParam() );
        int reuseInformationStartingFromTimeIndex = std::tr1::get<3>( GetParam() );

        assert( spaceMappingAlgorithm > -1 && spaceMappingAlgorithm < 3 );

        ASSERT_NEAR( tau, 0.01, 1.0e-13 );
        ASSERT_NEAR( kappa, 10, 1.0e-13 );
        ASSERT_TRUE( dx > 0 );

        // Initialize variables
        shared_ptr<TubeFlowFluidSolver> fluid;
        shared_ptr<TubeFlowSolidSolver> solid;
        shared_ptr<TubeFlowFluidSolver> fineModelFluid;
        shared_ptr<TubeFlowSolidSolver> fineModelSolid;
        shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        shared_ptr<PostProcessing> postProcessing;
        shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > solvers;
        shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

        solvers = shared_ptr< std::deque<shared_ptr<SpaceMappingSolver> > > ( new std::deque<shared_ptr<SpaceMappingSolver> >() );
        models = shared_ptr< std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > ( new std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );

        fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
        fineModelSolid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

        // Coarse model

        int maxUsedIterations = N;

        if ( parallel )
            maxUsedIterations *= 2;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

        shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
        shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
        shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
        shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, 5 ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 1.0e-3 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 1.0e-3 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        fluid.reset();
        solid.reset();
        multiLevelFluidSolver.reset();
        multiLevelSolidSolver.reset();
        convergenceMeasures.reset();
        multiLevelFsiSolver.reset();
        postProcessing.reset();

        // Fine model

        maxUsedIterations = 2 * N;

        if ( parallel )
            maxUsedIterations *= 2;

        fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, 2 * N, L, T, rho ) );
        solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, 2 * N ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 1 ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 1 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );
        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, 1.0e-2 * tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, 1.0e-2 * tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        shared_ptr<ImplicitMultiLevelFsiSolver> fineModel( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create manifold mapping object

        shared_ptr<SpaceMapping> spaceMapping;

        if ( spaceMappingAlgorithm == 0 )
            spaceMapping = shared_ptr<SpaceMapping>( new ManifoldMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian, true ) );

        if ( spaceMappingAlgorithm == 1 )
        {
            spaceMapping = shared_ptr<SpaceMapping>( new OutputSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );
        }

        if ( spaceMappingAlgorithm == 2 )
            spaceMapping = shared_ptr<SpaceMapping>( new AggressiveSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

        assert( spaceMapping );

        // Create manifold mapping solver
        shared_ptr<SpaceMappingSolver > spaceMappingSolver( new SpaceMappingSolver( fineModel, coarseModel, spaceMapping ) );

        solvers->push_back( spaceMappingSolver );
        models->push_back( coarseModel );
        models->push_back( fineModel );

        shared_ptr<SpaceMappingSolver> coarseModelSolver = spaceMappingSolver;

        fluid.reset();
        solid.reset();
        multiLevelFluidSolver.reset();
        multiLevelSolidSolver.reset();
        convergenceMeasures.reset();
        multiLevelFsiSolver.reset();
        postProcessing.reset();

        // Fine model 2

        maxUsedIterations = couplingGridSize;

        if ( parallel )
            maxUsedIterations *= 2;

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 2 ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToCouplingMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfFunction = shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
        rbfInterpolator = shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
        rbfInterpToMesh = shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, liveSelection, false, tol, tolLiveSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 2 ) );

        // Convergence measures
        convergenceMeasures = shared_ptr<std::list<shared_ptr<ConvergenceMeasure> > >( new std::list<shared_ptr<ConvergenceMeasure> > );

        convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, false, tol ) ) );

        if ( parallel )
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, false, tol ) ) );

        multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );
        postProcessing = shared_ptr<AndersonPostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        // Create space mapping object

        if ( spaceMappingAlgorithm == 0 )
            spaceMapping = shared_ptr<SpaceMapping>( new ManifoldMapping( fineModel, coarseModelSolver, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian, true ) );

        if ( spaceMappingAlgorithm == 1 )
            spaceMapping = shared_ptr<SpaceMapping>( new OutputSpaceMapping( fineModel, coarseModelSolver, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

        if ( spaceMappingAlgorithm == 2 )
            spaceMapping = shared_ptr<SpaceMapping>( new AggressiveSpaceMapping( fineModel, coarseModelSolver, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

        // Create manifold mapping solver
        spaceMappingSolver = shared_ptr<SpaceMappingSolver> ( new SpaceMappingSolver( fineModel, coarseModelSolver, spaceMapping ) );

        solvers->push_back( spaceMappingSolver );
        models->push_back( fineModel );

        // Create multi level manifold mapping solver

        solver = shared_ptr<fsi::MultiLevelSpaceMappingSolver> ( new fsi::MultiLevelSpaceMappingSolver( solvers, models, true ) );
    }

    virtual void TearDown()
    {
        solver.reset();
    }

    shared_ptr<fsi::MultiLevelSpaceMappingSolver> solver;
};

INSTANTIATE_TEST_CASE_P( testParameters, MultiLevelSpaceMappingSolverParametrizedTest, ::testing::Combine( Values( 2 ), Values( 30, 40 ), Values( 0, 1, 2 ), Values( 0, 5 ) ) );

TEST_P( MultiLevelSpaceMappingSolverParametrizedTest, run )
{
    for ( int i = 0; i < 5; i++ )
        solver->solveTimeStep();

    ASSERT_TRUE( solver->models->at( solver->models->size() - 1 )->fsi->allConverged );
    ASSERT_TRUE( solver->models->at( solver->models->size() - 1 )->fsi->fluid->isRunning() );
    ASSERT_EQ( solver->solvers->at( 0 )->spaceMapping->timeIndex, 5 );
    ASSERT_EQ( solver->solvers->at( 1 )->spaceMapping->timeIndex, 5 );
}
