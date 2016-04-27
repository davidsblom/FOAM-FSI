
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "SDC.H"
#include "Piston.H"
#include "Cos.H"
#include "Oscillator.H"
#include "gtest/gtest.h"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "ClenshawCurtis.H"
#include "UserDefinedNodes.H"

using namespace sdc;

TEST( UserDefinedNodes, run )
{
    std::vector<scalar> nodes ( 4, scalar( 0.0 ) );
    nodes[0] = 0;
    nodes[1] = 1./3.;
    nodes[2] = 2./3.;
    nodes[3] = 1.;

    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 5;

    endTime = 100;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-9;

    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature ( new fsi::quadrature::UserDefinedNodes<scalar>( nodes ) );

    std::shared_ptr<Piston> piston ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc ( new SDC( piston, quadrature, tol, 4, 10 ) );

    sdc->run();

    ASSERT_TRUE( sdc->isConverged() );
}

TEST( UserDefinedNodes, order )
{
    std::vector<scalar> nodes ( 4, scalar( 0.0 ) );
    nodes[0] = 0;
    nodes[1] = 1./3.;
    nodes[2] = 2./3.;
    nodes[3] = 1.;

    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 5;

    endTime = 100;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-9;

    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature1 ( new fsi::quadrature::UserDefinedNodes<scalar>( nodes ) );
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature2 ( new fsi::quadrature::UserDefinedNodes<scalar>( nodes ) );

    std::shared_ptr<Piston> piston1 ( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<Piston> piston2 ( new Piston( nbTimeSteps * 2, dt / 2., q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc1 ( new SDC( piston1, quadrature1, tol, 4, 10 ) );
    std::shared_ptr<SDC> sdc2 ( new SDC( piston2, quadrature1, tol, 4, 10 ) );

    sdc1->run();
    sdc2->run();

    fsi::vector solution1( 2 ), solution2( 2 ), f;
    piston1->getSolution( solution1, f );
    piston2->getSolution( solution2, f );

    scalar result1 = solution1( 1 );
    scalar result2 = solution2( 1 );

    scalar ref = piston1->referenceSolution( 100 );
    scalar error1 = std::abs( result1 - ref ) / std::abs( ref );
    scalar error2 = std::abs( result2 - ref ) / std::abs( ref );

    Info << "error1 = " << error1 << ", error2 = " << error2 << endl;

    scalar order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    ASSERT_GE( order, 4 );
}
