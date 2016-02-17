
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <mpi.h>
#include "gtest/gtest.h"
#include "dealiiSolidSolver.H"
#include "SDC.H"
#include "GaussRadau.H"

TEST( dealiiSolidSolver, sdc )
{
    using namespace dealiifsi;
    using namespace fsi;

    double time_step = 2.5e-3;
    double theta = 1;
    unsigned int degree = 1;
    unsigned int n_global_refines = 0;
    double gravity = 2;
    double distributed_load = 0;
    double rho = 1000;
    double final_time = 0.05;
    double nu = 0.4;
    double E = 1.4e6;

    std::shared_ptr<sdc::SDCSolver> solver( new dealiiSolidSolver<2> ( time_step, final_time, theta, degree, gravity, distributed_load, rho, E, nu, n_global_refines ) );

    int nbNodes = 2;
    scalar tol = 1.0e-15;
    int minSweeps = 1;
    int maxSweeps = 10;

    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

    sdc::SDC timeIntegrator( solver, quadrature, tol, minSweeps, maxSweeps );
    timeIntegrator.run();
}

TEST( dealiiSolidSolver, sdc_3 )
{
    using namespace dealiifsi;
    using namespace fsi;

    double time_step = 2.5e-3;
    double theta = 1;
    unsigned int degree = 1;
    unsigned int n_global_refines = 0;
    double gravity = 2;
    double distributed_load = 0;
    double rho = 1000;
    double final_time = 0.05;
    double nu = 0.4;
    double E = 1.4e6;

    std::shared_ptr<sdc::SDCSolver> solver( new dealiiSolidSolver<2> ( time_step, final_time, theta, degree, gravity, distributed_load, rho, E, nu, n_global_refines ) );

    int nbNodes = 3;
    scalar tol = 1.0e-10;
    int minSweeps = 1;
    int maxSweeps = 100;

    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

    sdc::SDC timeIntegrator( solver, quadrature, tol, minSweeps, maxSweeps );
    timeIntegrator.run();
}

TEST( dealiiSolidSolver, sdc_order )
{
    using namespace dealiifsi;
    using namespace fsi;

    double time_step = 2.5e-3;
    double theta = 1;
    unsigned int degree = 1;
    unsigned int n_global_refines = 0;
    double gravity = 2;
    double distributed_load = 0;
    double rho = 1000;
    double final_time = 0.05;
    double nu = 0.4;
    double E = 1.4e6;

    unsigned int nbComputations = 4;

    std::vector<double> solution_l2_norm( nbComputations );

    for ( unsigned int i = 0; i < nbComputations; ++i )
    {
        double dt = time_step / std::pow( 2, i );

        std::shared_ptr<sdc::SDCSolver> solver( new dealiiSolidSolver<2> ( dt, final_time, theta, degree, gravity, distributed_load, rho, E, nu, n_global_refines ) );

        int nbNodes = 2;
        scalar tol = 1.0e-10;
        int minSweeps = 1;
        int maxSweeps = 100;

        std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
        quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( nbNodes ) );

        sdc::SDC timeIntegrator( solver, quadrature, tol, minSweeps, maxSweeps );
        timeIntegrator.run();

        std::shared_ptr<dealiiSolidSolver<2> > elasticSolver = std::dynamic_pointer_cast<dealiiSolidSolver<2> > ( solver );
        fsi::vector solution( elasticSolver->getDOF() ), f( elasticSolver->getDOF() );
        elasticSolver->getSolution( solution, f );

        solution_l2_norm[i] = solution.norm();
    }

    std::vector<double> error( nbComputations - 1 );

    for ( unsigned int i = 0; i < solution_l2_norm.size(); ++i )
        std::cout << "l2norm = " << solution_l2_norm[i] << std::endl;

    for ( unsigned int i = 0; i < error.size(); ++i )
    {
        error[i] = std::abs( solution_l2_norm[i] - solution_l2_norm[nbComputations - 1] ) / std::abs( solution_l2_norm[nbComputations - 1] );

        std::cout << "error = " << error[i] << std::endl;
    }

    std::vector<double> order( nbComputations - 2 );

    for ( unsigned int i = 0; i < order.size(); ++i )
    {
        double dti = time_step / std::pow( 2, i );
        double dtinew = time_step / std::pow( 2, i + 1 );
        order[i] = std::log10( error[i + 1] ) - std::log10( error[i] );
        order[i] /= std::log10( dtinew ) - std::log10( dti );
        std::cout << "order = " << order[i] << std::endl;

        ASSERT_GE( order[i], 3 );
    }
}
