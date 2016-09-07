
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <chrono>

#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowSolidSolver.H"
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

using rbf::RBFFunctionInterface;
using rbf::RBFInterpolation;
using rbf::RBFCoarsening;
using rbf::TPSFunction;
using fsi::MultiLevelSolver;
using fsi::ConvergenceMeasure;
using fsi::RelativeConvergenceMeasure;
using fsi::MinIterationConvergenceMeasure;
using fsi::MultiLevelFsiSolver;
using fsi::PostProcessing;
using fsi::AndersonPostProcessing;
using fsi::ResidualRelativeConvergenceMeasure;
using fsi::AbsoluteConvergenceMeasure;

int main()
{
    int nbComputations = 6;
    int nbNodes = 6;
    std::vector<std::string> timeIntegrationSchemes = {
        "IDC", "SDIRK"
    };
    std::vector<std::string> sdirkSchemes = {
        "SDIRK2", "SDIRK3", "SDIRK4", "ESDIRK74PR", "SDIRK2PR"
    };

    for ( auto timeIntegrationSchemeString : timeIntegrationSchemes )
    {
        int nbSchemes = nbNodes;

        if ( timeIntegrationSchemeString == "SDIRK" )
            nbSchemes = sdirkSchemes.size();

        #pragma omp parallel for collapse(2), schedule(dynamic,1)

        for ( int iNodes = 0; iNodes < nbSchemes; iNodes++ )
        {
            for ( int iComputation = 0; iComputation < nbComputations; iComputation++ )
            {
                unsigned int nbNodes = iNodes + 1;
                unsigned int nbTimeSteps = std::pow( 2, iComputation );

                std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme;
                std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid;
                std::shared_ptr<MultiLevelFsiSolver> fsi;

                {
                    scalar r0 = 0.2;
                    scalar a0 = M_PI * r0 * r0;
                    scalar u0 = 0.1;
                    scalar p0 = 0;
                    scalar L = 1;
                    scalar T = 1;
                    scalar dt = T / nbTimeSteps;
                    scalar rho = 1.225;
                    scalar E = 490;
                    scalar h = 1.0e-3;
                    scalar cmk = std::sqrt( E * h / (2 * rho * r0) );

                    int N = 250;
                    bool parallel = false;
                    int extrapolation = 0;
                    scalar tol = 1.0e-5;
                    scalar absoluteTol = 1.0e-13;
                    int maxIter = 100;
                    scalar initialRelaxation = 1.0e-3;
                    int maxUsedIterations = 50;
                    int nbReuse = 0;

                    scalar singularityLimit = 1.0e-13;
                    int reuseInformationStartingFromTimeIndex = 0;
                    bool scaling = false;
                    bool updateJacobian = false;
                    scalar beta = 0.5;

                    fluid = std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> ( new tubeflow::SDCTubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
                    std::shared_ptr<tubeflow::SDCTubeFlowSolidSolver> solid( new tubeflow::SDCTubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

                    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;
                    std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
                    std::shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
                    std::shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

                    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new TPSFunction() );
                    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
                    rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

                    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new TPSFunction() );
                    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
                    rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

                    shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

                    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new TPSFunction() );
                    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
                    rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

                    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new TPSFunction() );
                    rbfInterpolator = std::shared_ptr<rbf::RBFInterpolation>( new RBFInterpolation( rbfFunction ) );
                    rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new RBFCoarsening( rbfInterpolator ) );

                    shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

                    std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
                    convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >( new std::list<std::shared_ptr<ConvergenceMeasure> > );

                    if ( timeIntegrationSchemeString == "IDC" )
                        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new ResidualRelativeConvergenceMeasure( 0, true, tol ) ) );

                    if ( timeIntegrationSchemeString == "SDIRK" )
                        convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new RelativeConvergenceMeasure( 0, true, absoluteTol ) ) );

                    convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure>( new AbsoluteConvergenceMeasure( 0, true, 0.1 * absoluteTol ) ) );

                    fsi = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

                    shared_ptr<PostProcessing> postProcessing( new AndersonPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

                    std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( fluid );
                    std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>( solid );

                    assert( sdcFluidSolver );
                    assert( sdcSolidSolver );

                    std::shared_ptr<fsi::SDCFsiSolver> fsiSolver( new fsi::SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing, extrapolation ) );

                    std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
                    quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::Uniform<scalar>( nbNodes ) );

                    if ( timeIntegrationSchemeString == "IDC" )
                        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( fsiSolver, quadrature, absoluteTol, 1, 50 ) );

                    if ( timeIntegrationSchemeString == "SDIRK" )
                    {
                        if ( iNodes + 1 > int( sdirkSchemes.size() ) )
                            continue;

                        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper( new sdc::AdaptiveTimeStepper( false ) );
                        std::string method = sdirkSchemes.at( iNodes );
                        timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::ESDIRK( fsiSolver, method, adaptiveTimeStepper ) );
                    }
                }

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
                ofstream data_u( label + "_data_u.log" );
                ofstream data_a( label + "_data_a.log" );
                ofstream data_p( label + "_data_p.log" );

                log_file << "label = " << label << std::endl;
                log_file << "nbNodes = " << nbNodes << std::endl;
                log_file << "nbTimeSteps = " << nbTimeSteps << std::endl;
                log_file << "timeIntegrationScheme = " << timeIntegrationSchemeString << std::endl;
                log_file << "nbIterations = " << fsi->nbIter << std::endl;
                log_file << "timing = " << elapsed_seconds.count() << std::endl;

                if ( timeIntegrationSchemeString == "SDIRK" )
                    log_file << "method = " << sdirkSchemes.at( iNodes ) << std::endl;

                data_u << std::setprecision( 20 ) << fluid->u << std::endl;
                data_a << std::setprecision( 20 ) << fluid->a << std::endl;
                data_p << std::setprecision( 20 ) << fluid->p << std::endl;

                log_file.close();
                data_a.close();
                data_u.close();
            }
        }
    }
}
