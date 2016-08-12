
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <chrono>

#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowBDFLinearizedSolidSolver.H"
#include "RBFCoarsening.H"
#include "SDC.H"
#include "ESDIRK.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "AbsoluteConvergenceMeasure.H"
#include "RelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "Uniform.H"
#include "AdaptiveTimeStepper.H"

int main()
{
    int nbComputations = 6;
    int nbNodes = 5;

    #pragma omp parallel for collapse(2), schedule(dynamic,1)

    for ( int iNodes = 0; iNodes < nbNodes; iNodes++ )
    {
        for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
        {
            unsigned int nbNodes = iNodes + 1;
            unsigned int nbTimeSteps = std::pow( 2, iComputation );

            std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme;
            std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid;
            std::shared_ptr<tubeflow::SDCTubeFlowBDFLinearizedSolidSolver> solid;
            std::shared_ptr<MultiLevelFsiSolver> fsi;

            {
                scalar r0 = 0.2;
                scalar a0 = M_PI * r0 * r0;
                scalar u0 = 0.1;
                scalar p0 = 0;
                scalar L = 1;
                scalar T = 1;
                scalar dt = T / nbTimeSteps;
                scalar rho_f = 1.225;
                scalar rho_s = 1.225;
                scalar E0 = 490;
                scalar G = 490;
                scalar h = 1.0e-3;
                scalar nu = 0.5;
                scalar cmk = std::sqrt( E0 * h / (2 * rho_f * r0) );

                int N = 250;
                bool parallel = false;
                int extrapolation = 0;
                int maxIter = 100;
                scalar initialRelaxation = 1.0e-3;
                int maxUsedIterations = 50;
                int nbReuse = 0;
                scalar tol = 1.0e-5;
                scalar absoluteTol = 1.0e-13;

                scalar singularityLimit = 1.0e-13;
                int reuseInformationStartingFromTimeIndex = 0;
                bool scaling = false;
                bool updateJacobian = false;
                scalar beta = 0.1;

                int timeOrder = 1;

                if ( nbNodes > 1 )
                    timeOrder = 2;

                fluid = std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> ( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
                solid = std::shared_ptr<tubeflow::SDCTubeFlowBDFLinearizedSolidSolver>( new tubeflow::SDCTubeFlowBDFLinearizedSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0, T, timeOrder ) );

                shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, 0, 0 ) );
                shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, 1, 0 ) );

                std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
                convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

                convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new ResidualRelativeConvergenceMeasure( 0, true, tol ) ) );
                convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new AbsoluteConvergenceMeasure( 0, true, 0.1 * absoluteTol ) ) );

                fsi = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

                shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

                std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( fluid );
                std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( solid );

                assert( sdcFluidSolver );
                assert( sdcSolidSolver );

                std::shared_ptr<fsi::SDCFsiSolver> fsiSolver( new fsi::SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing, extrapolation ) );

                std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::Uniform<scalar>( nbNodes ) );

                timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( fsiSolver, quadrature, absoluteTol, nbNodes, 50 ) );
            }

            assert( timeIntegrationScheme );
            assert( fluid );
            assert( solid );
            assert( fsi );

            std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
            start = std::chrono::high_resolution_clock::now();

            timeIntegrationScheme->run();

            end = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - start;

            std::string label = "IDC_BDF";
            label += "_nbNodes_" + std::to_string( nbNodes );
            label += "_nbTimeSteps_" + std::to_string( nbTimeSteps );

            ofstream log_file( label + ".log" );
            ofstream data_fluid_u( label + "_data_fluid_u.log" );
            ofstream data_fluid_a( label + "_data_fluid_a.log" );
            ofstream data_fluid_p( label + "_data_fluid_p.log" );
            ofstream data_solid_u( label + "_data_solid_u.log" );
            ofstream data_solid_r( label + "_data_solid_r.log" );

            log_file << "label = " << label << std::endl;
            log_file << "nbNodes = " << nbNodes << std::endl;
            log_file << "nbTimeSteps = " << nbTimeSteps << std::endl;
            log_file << "nbIterations = " << fsi->nbIter << std::endl;
            log_file << "timing = " << elapsed_seconds.count() << std::endl;

            data_fluid_u << std::setprecision( 20 ) << fluid->u << std::endl;
            data_fluid_a << std::setprecision( 20 ) << fluid->a << std::endl;
            data_fluid_p << std::setprecision( 20 ) << fluid->p << std::endl;
            data_solid_u << std::setprecision( 20 ) << solid->u << std::endl;
            data_solid_r << std::setprecision( 20 ) << solid->r << std::endl;

            log_file.close();
            data_fluid_u.close();
            data_fluid_a.close();
            data_fluid_p.close();
            data_solid_u.close();
            data_solid_r.close();
        }
    }
}
