
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ConvergenceMeasure.H"
#include "FsiSolver.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MinIterationConvergenceMeasure.H"
#include "MonolithicFsiSolver.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"
#include "AndersonPostProcessing.H"
#include <unsupported/Eigen/NumericalDiff>

using namespace tubeflow;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

class AndersonPostProcessingParametrizedTest : public TestWithParam< std::tr1::tuple<bool, int, int, int, bool> >
{
protected:

  virtual void SetUp()
  {
    // Physical settings
    double r0 = 0.2;
    double a0 = M_PI * r0 * r0;
    double u0 = 0.1;
    double p0 = 0;
    double dt = 0.1;
    int N = 5;
    double L = 1;
    double T = 10;
    double dx = L / N;
    double rho = 1.225;
    double E = 490;
    double h = 1.0e-3;
    double cmk = std::sqrt( E * h / (2 * rho * r0) );
    double c0 = std::sqrt( cmk * cmk - p0 / (2 * rho) );
    double kappa = c0 / u0;
    double tau = u0 * dt / L;

    // Computational settings
    double tol = 1.0e-7;
    int maxIter = 50;
    double initialRelaxation = 1.0e-3;
    double singularityLimit = 1.0e-11;
    int reuseInformationStartingFromTimeIndex = 0;
    bool scaling = false;
    double beta = 1;

    // Parametrized settings
    bool parallel = std::tr1::get<0>( GetParam() );
    int nbReuse = std::tr1::get<1>( GetParam() );
    int extrapolation = std::tr1::get<2>( GetParam() );
    int minIter = std::tr1::get<3>( GetParam() );
    bool updateJacobian = std::tr1::get<4>( GetParam() );

    int maxUsedIterations = N;

    if ( parallel )
      maxUsedIterations *= 2;

    ASSERT_NEAR( tau, 0.01, 1.0e-13 );
    ASSERT_NEAR( kappa, 10, 1.0e-13 );
    ASSERT_TRUE( dx > 0 );

    shared_ptr<TubeFlowFluidSolver> fluid( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    shared_ptr<TubeFlowSolidSolver> solid( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

    shared_ptr<RBFFunctionInterface> rbfFunction;
    shared_ptr<RBFInterpolation> rbfInterpolator;
    shared_ptr<RBFCoarsening> rbfInterpToCouplingMesh;
    shared_ptr<RBFCoarsening> rbfInterpToMesh;




    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

    shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

    rbfFunction = shared_ptr<RBFFunctionInterface>( new TPSFunction() );
    rbfInterpolator = shared_ptr<RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = shared_ptr<RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

    shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

    // Convergence measures
    std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

    convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, minIter ) ) );
    convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, tol ) ) );

    if ( parallel )
      convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 1, tol ) ) );

    shared_ptr<MultiLevelFsiSolver> fsi( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );
    shared_ptr<AndersonPostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );
    solver = new ImplicitMultiLevelFsiSolver( fsi, postProcessing );
    monolithicSolver = new MonolithicFsiSolver( a0, u0, p0, dt, cmk, N, L, T, rho );
  }

  virtual void TearDown()
  {
    delete solver;
    delete monolithicSolver;
  }

  ImplicitMultiLevelFsiSolver * solver;
  MonolithicFsiSolver * monolithicSolver;
};

INSTANTIATE_TEST_CASE_P( testParameters, AndersonPostProcessingParametrizedTest, ::testing::Combine( Bool(), Values( 0, 1, 4 ), Values( 0, 1, 2 ), Values( 1 ), Bool() ) );

TEST_P( AndersonPostProcessingParametrizedTest, object )
{
  ASSERT_TRUE( true );
}

TEST_P( AndersonPostProcessingParametrizedTest, run )
{
  solver->run();
  monolithicSolver->run();

  double tol = 1.0e-5;
  ASSERT_TRUE( solver->fsi->allConverged );
  ASSERT_NEAR( solver->fsi->fluid->data.norm(), monolithicSolver->pn.norm(), tol );
  ASSERT_NEAR( solver->fsi->solid->data.norm(), monolithicSolver->an.norm(), tol );
  ASSERT_TRUE( monolithicSolver->an.norm() > 0 );
  ASSERT_TRUE( monolithicSolver->pn.norm() > 0 );
}
