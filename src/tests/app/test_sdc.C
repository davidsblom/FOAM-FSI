
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "gtest/gtest.h"

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

class Piston
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

  void solve();

  void solveTimeStep(
    const double t,
    const double dt,
    const Eigen::VectorXd & qold,
    const Eigen::VectorXd rhs,
    Eigen::VectorXd & f,
    Eigen::VectorXd & result
    );

  void evaluate(
    const Eigen::VectorXd & q,
    const double t,
    Eigen::VectorXd & f
    );

  int nbTimeSteps;
  double dt;
  double q0;
  double c1;
  double c2;
  double As;
  double Ac;
  double omega;
  int N;
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
  nbTimeSteps( nbTimeSteps ),
  dt( dt ),
  q0( q0 ),
  c1( q0 + Ac / std::pow( omega, 2 ) ),
  c2( qdot0 + As / omega ),
  As( As ),
  Ac( Ac ),
  omega( omega ),
  N( 2 )
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

void Piston::evaluate(
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

void Piston::solveTimeStep(
  const double t,
  const double dt,
  const Eigen::VectorXd & qold,
  const Eigen::VectorXd rhs,
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
}

class PistonTest : public::testing::Test
{
protected:

  virtual void SetUp()
  {
    int nbTimeSteps;
    double dt, q0, qdot0, As, Ac, omega, endTime;

    nbTimeSteps = 1000;
    endTime = 100;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;

    piston = new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega );
  }

  virtual void TearDown()
  {
    delete piston;
    piston = NULL;
  }

  Piston * piston;
};

TEST_F( PistonTest, object )
{
  ASSERT_TRUE( true );
}

TEST_F( PistonTest, referenceSolution )
{
  double result = piston->referenceSolution( 100 );

  ASSERT_NEAR( result, -35.5953231178, 1.0e-11 );
}
