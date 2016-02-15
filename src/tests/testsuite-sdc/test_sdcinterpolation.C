
#include <Eigen/Dense>
#include "SDC.H"
#include "gtest/gtest.h"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Piston.H"

TEST( SDC, interpolation )
{
    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 1;

    endTime = 10;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-17;

    int nbNodes = 10;
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

    std::shared_ptr<Piston> piston( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc( new SDC( piston, quadrature, tol, 5, 10 * nbNodes ) );

    sdc->run();

    fsi::vector solution( 2 ), f;
    piston->getSolution( solution, f );
    scalar result = solution( 1 );
    scalar ref = piston->referenceSolution( endTime );
    scalar error = std::abs( result - ref ) / std::abs( ref );

    std::vector<scalar> nodes = {
        0.856, 0.365, 1.0
    };

    const fsi::matrix interp = sdc->data->integrate( nodes, dt );

    for ( unsigned int i = 0; i < nodes.size(); i++ )
    {
        scalar ref = piston->referenceSolution( nodes[i] * endTime );
        scalar errori = abs( interp( i, 1 ) - ref ) / abs( ref );
        std::cout << "error" << i << " = " << errori << std::endl;

        if ( i == 2 )
            ASSERT_NEAR( error, errori, 1.0e-13 );
    }
}

TEST( SDC, interpolationOrder )
{
    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 1;

    endTime = 1;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-17;

    int nbNodes = 2;
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

    std::shared_ptr<Piston> piston1( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc1( new SDC( piston1, quadrature, tol, 5, 10 * nbNodes ) );
    std::shared_ptr<Piston> piston2( new Piston( nbTimeSteps * 2, dt / 2.0, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc2( new SDC( piston2, quadrature, tol, 5, 10 * nbNodes ) );

    sdc1->run();
    sdc2->run();

    std::vector<scalar> nodes1 = {
        0.8
    };
    std::vector<scalar> nodes2 = {
        0.6
    };

    const fsi::matrix interp1 = sdc1->data->integrate( nodes1, dt );
    const fsi::matrix interp2 = sdc2->data->integrate( nodes2, 0.5 * dt );

    scalar ref1 = piston1->referenceSolution( nodes1[0] * endTime );
    scalar ref2 = piston2->referenceSolution( nodes1[0] * endTime );

    scalar error1 = std::abs( interp1( 0, 1 ) - ref1 ) / std::abs( ref1 );
    scalar error2 = std::abs( interp2( 0, 1 ) - ref2 ) / std::abs( ref2 );

    scalar order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    ASSERT_GE( order, 3 );
}

TEST( SDC, interpolationGaussLobatto )
{
    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 1;

    endTime = 1;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-17;

    int nbNodes = 2;
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussLobatto<scalar>( nbNodes ) );

    std::shared_ptr<Piston> piston1( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc1( new SDC( piston1, quadrature, tol, 5, 10 * nbNodes ) );
    std::shared_ptr<Piston> piston2( new Piston( nbTimeSteps * 2, dt / 2.0, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc2( new SDC( piston2, quadrature, tol, 5, 10 * nbNodes ) );

    sdc1->run();
    sdc2->run();

    std::vector<scalar> nodes1 = {
        0.8
    };
    std::vector<scalar> nodes2 = {
        0.6
    };

    const fsi::matrix interp1 = sdc1->data->integrate( nodes1, dt );
    const fsi::matrix interp2 = sdc2->data->integrate( nodes2, 0.5 * dt );

    scalar ref1 = piston1->referenceSolution( nodes1[0] * endTime );
    scalar ref2 = piston2->referenceSolution( nodes1[0] * endTime );

    scalar error1 = std::abs( interp1( 0, 1 ) - ref1 ) / std::abs( ref1 );
    scalar error2 = std::abs( interp2( 0, 1 ) - ref2 ) / std::abs( ref2 );

    scalar order = ( std::log10( error1 ) - std::log10( error2 ) ) / ( std::log10( nbTimeSteps * 2 ) - std::log10( nbTimeSteps ) );

    ASSERT_GE( order, 4 );
}

TEST( SDC, polynomialInterpolationGaussLobatto )
{
    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 1;

    endTime = 10;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-17;

    int nbNodes = 10;
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussLobatto<scalar>( nbNodes ) );

    std::shared_ptr<Piston> piston( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc( new SDC( piston, quadrature, tol, 5, 10 * nbNodes ) );

    fsi::matrix functions( quadrature->get_num_nodes(), 1 );
    int i = 0;

    for ( auto node : quadrature->get_nodes() )
    {
        functions( i, 0 ) = std::sin( node );
        i++;
    }

    std::vector<scalar> nodes = {
        0, 0.5, 0.3645, quadrature->get_nodes()[1]
    };

    const fsi::matrix interp = sdc->data->interpolate( functions, nodes );

    i = 0;

    for ( auto node : nodes )
    {
        scalar ref = std::sin( node );
        scalar error = abs( interp( i, 0 ) - ref ) / abs( ref + SMALL );
        ASSERT_LE( error, 1.0e-10 );
        i++;
    }
}

TEST( SDC, polynomialInterpolationGaussRadau )
{
    scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

    int nbTimeSteps = 1;

    endTime = 10;
    dt = endTime / nbTimeSteps;
    As = 100;
    Ac = As;
    omega = 1;
    q0 = -As;
    qdot0 = -As;
    tol = 1.0e-17;

    int nbNodes = 10;
    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

    std::shared_ptr<Piston> piston( new Piston( nbTimeSteps, dt, q0, qdot0, As, Ac, omega ) );
    std::shared_ptr<SDC> sdc( new SDC( piston, quadrature, tol, 5, 10 * nbNodes ) );

    fsi::matrix functions( quadrature->get_num_nodes(), 1 );
    int i = 0;

    for ( auto node : quadrature->get_nodes() )
    {
        functions( i, 0 ) = std::sin( node );
        i++;
    }

    std::vector<scalar> nodes = {
        0, 0.5, 0.3645, quadrature->get_nodes()[1]
    };

    const fsi::matrix interp = sdc->data->interpolate( functions, nodes );

    i = 0;

    for ( auto node : nodes )
    {
        scalar ref = std::sin( node );
        scalar error = abs( interp( i, 0 ) - ref ) / abs( ref + SMALL );
        ASSERT_LE( error, 1.0e-10 );
        i++;
    }
}
