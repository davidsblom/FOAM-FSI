
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <chrono>

#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowLinearSolidSolver.H"
#include "RBFCoarsening.H"
#include "SDC.H"
#include "ESDIRK.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "RelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "Uniform.H"
#include "AdaptiveTimeStepper.H"

int main()
{
    int nbComputations = 9;
    int nbNodes = 6;
    std::vector<std::string> timeIntegrationSchemes = {
        "IDC", "SDIRK"
    };
    std::vector<std::string> sdirkSchemes = {
        "SDIRK2", "SDIRK3", "SDIRK4", "ESDIRK3", "ESDIRK4", "ESDIRK5", "ESDIRK53PR", "ESDIRK63PR", "ESDIRK74PR"
    };

    for ( auto timeIntegrationSchemeString : timeIntegrationSchemes )
    {
        int nbSchemes = nbNodes;

        if ( timeIntegrationSchemeString == "SDIRK" )
            nbSchemes = sdirkSchemes.size();

        for ( int iNodes = 0; iNodes < nbSchemes; iNodes++ )
        {
            for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
            {
                unsigned int nbNodes = iNodes + 1;
                unsigned int nbTimeSteps = std::pow( 2, iComputation );

                std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme;
                std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid;
                std::shared_ptr<tubeflow::SDCTubeFlowLinearSolidSolver> solid;
                std::shared_ptr<MultiLevelFsiSolver> fsi;

                {
                    scalar r0 = 3.0e-3;
                    scalar h = 3.0e-4;
                    scalar L = 0.126;
                    scalar rho_s = 1000;
                    scalar E0 = 4.0e5;
                    scalar G = 4.0e5;
                    scalar nu = 0.5;

                    scalar a0 = M_PI * r0 * r0;
                    scalar u0 = 0.26;
                    scalar p0 = 0;
                    int N = 100;
                    scalar T = 1;
                    scalar dt = T / nbTimeSteps;
                    scalar rho_f = 1060;
                    scalar E = 490;
                    scalar cmk = std::sqrt( E * h / (2 * rho_f * r0) );

                    bool parallel = false;
                    int extrapolation = 0;
                    scalar tol = 1.0e-3;
                    int maxIter = 20;
                    scalar initialRelaxation = 1.0e-3;
                    int maxUsedIterations = 50;
                    int nbReuse = 0;

                    scalar singularityLimit = 1.0e-13;
                    int reuseInformationStartingFromTimeIndex = 0;
                    bool scaling = false;
                    bool updateJacobian = false;
                    scalar beta = 0.5;
                    int minIter = 5;

                    fluid = std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> ( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho_f ) );
                    solid = std::shared_ptr<tubeflow::SDCTubeFlowLinearSolidSolver>( new tubeflow::SDCTubeFlowLinearSolidSolver( N, nu, rho_s, h, L, dt, G, E0, r0 ) );

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

                    std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
                    convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

                    if ( timeIntegrationSchemeString == "IDC" )
                        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new ResidualRelativeConvergenceMeasure( 0, false, tol ) ) );

                    if ( timeIntegrationSchemeString == "SDIRK" )
                        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, false, tol ) ) );

                    convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new MinIterationConvergenceMeasure( 0, false, minIter ) ) );

                    fsi = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

                    shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

                    std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( fluid );
                    std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( solid );

                    assert( sdcFluidSolver );
                    assert( sdcSolidSolver );

                    std::shared_ptr<fsi::SDCFsiSolver> fsiSolver( new fsi::SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing ) );

                    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
                    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::Uniform<scalar>( nbNodes ) );

                    if ( timeIntegrationSchemeString == "IDC" )
                        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( fsiSolver, quadrature, 1.0e-13, 1, 50 ) );

                    if ( timeIntegrationSchemeString == "SDIRK" )
                    {
                        if ( iNodes + 1 > int( sdirkSchemes.size() ) )
                            continue;

                        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper( new sdc::AdaptiveTimeStepper( false ) );
                        std::string method = sdirkSchemes.at( iNodes );
                        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::ESDIRK( fsiSolver, method, adaptiveTimeStepper ) );
                    }
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

                std::string label = timeIntegrationSchemeString;

                if ( timeIntegrationSchemeString == "SDIRK" )
                    label += "_" + sdirkSchemes.at( iNodes );

                if ( timeIntegrationSchemeString == "IDC" )
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
                log_file << "timeIntegrationScheme = " << timeIntegrationSchemeString << std::endl;
                log_file << "nbIterations = " << fsi->nbIter << std::endl;
                log_file << "timing = " << elapsed_seconds.count() << std::endl;

                if ( timeIntegrationSchemeString == "SDIRK" )
                    log_file << "method = " << sdirkSchemes.at( iNodes ) << std::endl;

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
}
