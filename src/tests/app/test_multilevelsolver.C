
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MultiLevelSolver.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"
#include "gtest/gtest.h"

using namespace tubeflow;

class MultiLevelSolverTest : public::testing::Test
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
    int couplingGridSize = 2 * N;
    double L = 1;
    double T = 10;
    double rho = 1.225;
    double E = 490;
    double h = 1.0e-3;
    double cmk = std::sqrt( E * h / (2 * rho * r0) );

    // Initialize variables
    fluid = shared_ptr<TubeFlowFluidSolver>( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
    fluidFine = shared_ptr<TubeFlowFluidSolver>( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, couplingGridSize, L, T, rho ) );
    solid = shared_ptr<TubeFlowSolidSolver>( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
    solidFine = shared_ptr<TubeFlowSolidSolver>( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, couplingGridSize ) );

    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
    std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
    std::shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
    std::shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;
    bool coarsening = false;
    double tol = 1.0e-6;
    int coarseningMinPoints = 200;
    int coarseningMaxPoints = 2000;

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    fluidSolver = new MultiLevelSolver( fluid, fluidFine, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    fluidFineSolver = new MultiLevelSolver( fluidFine, fluidFine, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    solidSolver = new MultiLevelSolver( solid, solidFine, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );
    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction ) );
    rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, tol, coarseningMinPoints, coarseningMaxPoints ) );

    solidFineSolver = new MultiLevelSolver( solidFine, solidFine, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 );
  }

  virtual void TearDown()
  {
    fluid.reset();
    fluidFine.reset();
    solid.reset();
    solidFine.reset();
    delete fluidSolver;
    delete fluidFineSolver;
    delete solidSolver;
    delete solidFineSolver;
  }

  shared_ptr<TubeFlowFluidSolver> fluid;
  shared_ptr<TubeFlowFluidSolver> fluidFine;
  shared_ptr<TubeFlowSolidSolver> solid;
  shared_ptr<TubeFlowSolidSolver> solidFine;
  MultiLevelSolver * fluidSolver;
  MultiLevelSolver * fluidFineSolver;
  MultiLevelSolver * solidSolver;
  MultiLevelSolver * solidFineSolver;
};

TEST_F( MultiLevelSolverTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( MultiLevelSolverTest, fluid )
{
  int N = 5;
  matrix input( N, 1 ), output( N, 1 );
  input.fill( fluid->a0 );
  output.setZero();

  fluid->initTimeStep();
  fluid->solve( input, output );
}

TEST_F( MultiLevelSolverTest, fluidSolver )
{
  int N = 10;
  matrix input( N, 1 ), output( N, 1 );
  input.fill( fluid->a0 );
  output.setZero();

  fluid->initTimeStep();
  fluidSolver->updateInterpolationMatrices();
  fluidSolver->solve( input, output );
}

TEST_F( MultiLevelSolverTest, fluidFineSolver )
{
  int N = 10;
  matrix input( N, 1 ), output( N, 1 );
  input.fill( fluid->a0 );
  output.setZero();

  fluidFine->initTimeStep();
  fluidFineSolver->updateInterpolationMatrices();
  fluidFineSolver->solve( input, output );
}

TEST_F( MultiLevelSolverTest, solid )
{
  int N = 5;
  matrix input( N, 1 ), output( N, 1 );
  input.fill( fluid->p0 );
  output.setZero();

  solid->initTimeStep();
  solid->solve( input, output );
}

TEST_F( MultiLevelSolverTest, solidSolver )
{
  int N = 10;
  matrix input( N, 1 ), output( N, 1 );
  input.fill( fluid->p0 );
  output.setZero();

  solid->initTimeStep();
  solidSolver->updateInterpolationMatrices();
  solidSolver->solve( input, output );
}

TEST_F( MultiLevelSolverTest, solidFineSolver )
{
  int N = 10;
  matrix input( N, 1 ), output( N, 1 );
  input.fill( fluid->p0 );
  output.setZero();

  solidFine->initTimeStep();
  solidFineSolver->updateInterpolationMatrices();
  solidFineSolver->solve( input, output );
}

TEST_F( MultiLevelSolverTest, gridFluid )
{
  fluid->calcGrid();
  double tol = 1.0e-14;
  ASSERT_NEAR( fluid->grid( 0, 0 ), 0.1, tol );
  ASSERT_NEAR( fluid->grid( 1, 0 ), 0.3, tol );
  ASSERT_NEAR( fluid->grid( 2, 0 ), 0.5, tol );
  ASSERT_NEAR( fluid->grid( 3, 0 ), 0.7, tol );
  ASSERT_NEAR( fluid->grid( 4, 0 ), 0.9, tol );
  ASSERT_EQ( fluid->grid.cols(), 1 );
  ASSERT_EQ( fluid->grid.rows(), 5 );
}

TEST_F( MultiLevelSolverTest, gridSolid )
{
  solid->calcGrid();
  double tol = 1.0e-14;
  ASSERT_NEAR( solid->grid( 0, 0 ), 0.1, tol );
  ASSERT_NEAR( solid->grid( 1, 0 ), 0.3, tol );
  ASSERT_NEAR( solid->grid( 2, 0 ), 0.5, tol );
  ASSERT_NEAR( solid->grid( 3, 0 ), 0.7, tol );
  ASSERT_NEAR( solid->grid( 4, 0 ), 0.9, tol );
  ASSERT_EQ( solid->grid.cols(), 1 );
  ASSERT_EQ( solid->grid.rows(), 5 );
}

TEST_F( MultiLevelSolverTest, fineGridFluid )
{
  fluidFine->calcGrid();
  double tol = 1.0e-14;
  ASSERT_NEAR( fluidFine->grid( 0, 0 ), 0.05, tol );
  ASSERT_NEAR( fluidFine->grid( 1, 0 ), 0.15, tol );
  ASSERT_NEAR( fluidFine->grid( 2, 0 ), 0.25, tol );
  ASSERT_NEAR( fluidFine->grid( 3, 0 ), 0.35, tol );
  ASSERT_NEAR( fluidFine->grid( 4, 0 ), 0.45, tol );
  ASSERT_NEAR( fluidFine->grid( 5, 0 ), 0.55, tol );
  ASSERT_NEAR( fluidFine->grid( 6, 0 ), 0.65, tol );
  ASSERT_NEAR( fluidFine->grid( 7, 0 ), 0.75, tol );
  ASSERT_NEAR( fluidFine->grid( 8, 0 ), 0.85, tol );
  ASSERT_NEAR( fluidFine->grid( 9, 0 ), 0.95, tol );
  ASSERT_EQ( fluidFine->grid.cols(), 1 );
  ASSERT_EQ( fluidFine->grid.rows(), 10 );
}

TEST_F( MultiLevelSolverTest, fineGridSolid )
{
  solidFine->calcGrid();
  double tol = 1.0e-14;
  ASSERT_NEAR( solidFine->grid( 0, 0 ), 0.05, tol );
  ASSERT_NEAR( solidFine->grid( 1, 0 ), 0.15, tol );
  ASSERT_NEAR( solidFine->grid( 2, 0 ), 0.25, tol );
  ASSERT_NEAR( solidFine->grid( 3, 0 ), 0.35, tol );
  ASSERT_NEAR( solidFine->grid( 4, 0 ), 0.45, tol );
  ASSERT_NEAR( solidFine->grid( 5, 0 ), 0.55, tol );
  ASSERT_NEAR( solidFine->grid( 6, 0 ), 0.65, tol );
  ASSERT_NEAR( solidFine->grid( 7, 0 ), 0.75, tol );
  ASSERT_NEAR( solidFine->grid( 8, 0 ), 0.85, tol );
  ASSERT_NEAR( solidFine->grid( 9, 0 ), 0.95, tol );
  ASSERT_EQ( solidFine->grid.cols(), 1 );
  ASSERT_EQ( solidFine->grid.rows(), 10 );
}
