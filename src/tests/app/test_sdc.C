
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "SDC.H"
#include "gtest/gtest.h"

using namespace sdc;
using::testing::TestWithParam;
using::testing::Bool;
using::testing::Values;
using::testing::Combine;

/*
 * Implements the piston problem for which a exact solution exists.
 * With the exact solution, the convergence orders can be confirmed
 * exactly for this test case.
 *
 * Reference: A. H. van Zuijlen, Efficient Higher Order Time
 * Integration of Partitioned Systems, PhD thesis, page 22.
 *
 * Notes Gaussian quadrature:
 * http://undergraduate.csse.uwa.edu.au/units/GENG2140/Mech2012/Gaussian_Quadrature.pdf
 * http://en.wikipedia.org/wiki/Gaussian_quadrature
 * http://en.wikipedia.org/wiki/Legendre_polynomials
 */

class Piston : public SDCSolver
{
public:

  Piston(
    int nbTimeSteps,
    double dt,
    double q0,
    double qdot0,
    double As,
    double Ac,
    double omega
    );

  ~Piston();

  double referenceSolution( double t );

  void run();

  virtual void finalizeTimeStep(){}

  virtual void getSolution( Eigen::VectorXd & solution );

  virtual void initTimeStep(){}

  virtual bool isRunning();

  virtual void implicitSolve(
    bool corrector,
    const int k,
    const double t,
    const double dt,
    const Eigen::VectorXd & qold,
    const Eigen::VectorXd & rhs,
    Eigen::VectorXd & f,
    Eigen::VectorXd & result
    );

  virtual int getDOF();

  virtual double getTimeStep();

  virtual void evaluateFunction(
    const int k,
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
    );

  virtual void setDeltaT( double dt ){}

  virtual void setNumberOfStages( int k ){}

  virtual void nextTimeStep(){}

  virtual int getNbCells()
  {
    return 2;
  }

  virtual double getScalingFactor()
  {
    return 1;
  }

  int nbTimeSteps;
  double dt;
  double q0;
  double c1;
  double c2;
  double As;
  double Ac;
  double omega;
  int N;
  double q;
  double qdot;
  double t;
};

Piston::Piston(
  int nbTimeSteps,
  double dt,
  double q0,
  double qdot0,
  double As,
  double Ac,
  double omega
  )
  :
  SDCSolver(),
  nbTimeSteps( nbTimeSteps ),
  dt( dt ),
  q0( q0 ),
  c1( q0 + Ac / std::pow( omega, 2 ) ),
  c2( qdot0 + As / omega ),
  As( As ),
  Ac( Ac ),
  omega( omega ),
  N( 2 ),
  q( q0 ),
  qdot( qdot0 ),
  t( 0 )
{
  assert( nbTimeSteps > 0 );
  assert( dt > 0 );
}

Piston::~Piston()
{}

double Piston::referenceSolution( double t )
{
  double result = c1;
  result += c2 * t;
  result += -As / std::pow( omega, 2 ) * std::sin( omega * t );
  result += -Ac / std::pow( omega, 2 ) * std::cos( omega * t );

  return result;
}

void Piston::evaluateFunction(
  const int k,
  const Eigen::VectorXd & q,
  const double t,
  Eigen::VectorXd & f
  )
{
  assert( f.rows() == 2 );

  f( 0 ) = As * std::sin( omega * t );
  f( 0 ) += Ac * std::cos( omega * t );
  f( 1 ) = q( 0 );
}

int Piston::getDOF()
{
  return 2;
}

void Piston::getSolution( Eigen::VectorXd & solution )
{
  assert( solution.rows() == 2 );
  solution( 0 ) = qdot;
  solution( 1 ) = q;
}

double Piston::getTimeStep()
{
  return dt;
}

bool Piston::isRunning()
{
  if ( t > nbTimeSteps * dt )
    return false;

  return std::abs( t - nbTimeSteps * dt ) > 1.0e-10;
}

void Piston::run()
{
  Eigen::VectorXd q( nbTimeSteps + 1 ), qdot( nbTimeSteps + 1 ), qold( 2 ), f( 2 ), rhs( 2 ), result( 2 );

  q( 0 ) = -Ac;
  qdot( 0 ) = q( 0 );

  for ( int i = 1; i < nbTimeSteps + 1; i++ )
  {
    double t = dt * i;

    qold << qdot( i - 1 ), q( i - 1 );

    f.setZero();
    rhs.setZero();

    implicitSolve( false, 0, t, dt, qold, rhs, f, result );

    qdot( i ) = result( 0 );
    q( i ) = result( 1 );
  }
}

void Piston::implicitSolve(
  bool corrector,
  const int k,
  const double t,
  const double dt,
  const Eigen::VectorXd & qold,
  const Eigen::VectorXd & rhs,
  Eigen::VectorXd & f,
  Eigen::VectorXd & result
  )
{
  assert( f.rows() == 2 );
  assert( rhs.rows() == 2 );
  assert( result.rows() == 2 );

  f( 0 ) = As * std::sin( omega * t );
  f( 0 ) += Ac * std::cos( omega * t );

  // qdot
  result( 0 ) = qold( 0 ) + dt * f( 0 ) + rhs( 0 );

  // q
  result( 1 ) = qold( 1 ) + std::pow( dt, 2 ) * f( 0 ) + dt * qold( 0 ) + dt * rhs( 0 ) + rhs( 1 );

  f( 1 ) = result( 0 );

  qdot = result( 0 );
  q = result( 1 );

  this->t = t;
}

class SDCTest : public TestWithParam< std::tr1::tuple<int, int> >
{
protected:

  virtual void SetUp()
  {
    double dt, q0, qdot0, As, Ac, omega, endTime;

    int nbNodes = std::tr1::get<0>( GetParam() );
    int nbTimeSteps = std::tr1::get<1>( GetParam() );

    endTime = 100;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;

    piston = std::shared_ptr<Piston> ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    sdc = std::shared_ptr<SDC> ( new SDC( piston, nbNodes, 1.0e-13 ) );
  }

  virtual void TearDown()
  {
    piston.reset();
    sdc.reset();
  }

  std::shared_ptr<Piston> piston;
  std::shared_ptr<SDC> sdc;
};

INSTANTIATE_TEST_CASE_P( testParameters, SDCTest, ::testing::Combine( Values( 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 ), Values( 10, 20, 40, 80, 100, 160 ) ) );

TEST_P( SDCTest, object )
{
  ASSERT_TRUE( true );
}

TEST_P( SDCTest, referenceSolution )
{
  double result = piston->referenceSolution( 100 );

  ASSERT_NEAR( result, -35.5953231178, 1.0e-11 );
}

TEST_P( SDCTest, solve )
{
  piston->run();

  Eigen::VectorXd solution( 2 );
  piston->getSolution( solution );

  double result = solution( 1 );

  int nbTimeSteps = std::tr1::get<1>( GetParam() );

  if ( nbTimeSteps == 10 )
    ASSERT_NEAR( result, -75814.5607609, 1.0e-8 );

  ASSERT_NEAR( piston->t, 100, 1.0e-10 );
}

TEST_P( SDCTest, evaluateFunction )
{
  Eigen::VectorXd q( 2 ), f( 2 );
  double t;

  q << -100, -100;
  t = 1;

  piston->evaluateFunction( 0, q, t, f );

  ASSERT_NEAR( f( 0 ), 138.17732907, 1.0e-8 );
  ASSERT_NEAR( f( 1 ), -100, 1.0e-9 );
}

TEST_P( SDCTest, solveTimeStep )
{
  sdc->solveTimeStep( 0 );
}

TEST_P( SDCTest, run )
{
  sdc->run();

  Eigen::VectorXd solution( 2 );
  piston->getSolution( solution );

  double result = solution( 1 );
  double ref = piston->referenceSolution( 100 );
  double error = std::abs( result - ref ) / std::abs( ref );

  int nbTimeSteps = std::tr1::get<1>( GetParam() );
  int nbNodes = std::tr1::get<0>( GetParam() );

  ASSERT_NEAR( piston->t, 100, 1.0e-10 );

  if ( nbTimeSteps == 100 && nbNodes > 4 )
    ASSERT_NEAR( error, 1.0e-5, 1.0e-4 );

  if ( nbTimeSteps == 10 && nbNodes == 3 )
    ASSERT_NEAR( error, 1.72134678e+03, 1.0e-5 );

  if ( nbTimeSteps == 20 && nbNodes == 3 )
    ASSERT_NEAR( error, 5.38712906e+02, 1.0e-6 );

  if ( nbTimeSteps == 40 && nbNodes == 3 )
    ASSERT_NEAR( error, 2.62943671e+01, 1.0e-7 );

  if ( nbTimeSteps == 80 && nbNodes == 3 )
    ASSERT_NEAR( error, 2.74294224e+00, 1.0e-8 );

  if ( nbTimeSteps == 160 && nbNodes == 3 )
    ASSERT_NEAR( error, 3.25178345e-01, 1.0e-9 );
}
