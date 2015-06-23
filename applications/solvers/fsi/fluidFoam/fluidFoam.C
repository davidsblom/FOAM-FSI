
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>
#include <yaml-cpp/yaml.h>

#include "FluidSolver.H"
#include "CoupledFluidSolver.H"
#include "CompressibleFluidSolver.H"
#include "SDCFluidSolver.H"
#include "SDC.H"
#include "SDCLaplacianSolver.H"
#include "ESDIRK.H"
#include "AdaptiveTimeStepper.H"

int main(
    int argc,
    char * argv[]
    )
{
    std::shared_ptr<argList> args( new argList( argc, argv ) );

    if ( !args->checkRootCase() )
    {
        FatalError.exit();
    }

    std::shared_ptr<Time> runTime( new Time
        (
            Time::controlDictName,
            args->rootPath(),
            args->caseName()
        ) );

    string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";

    YAML::Node config = YAML::LoadFile( filename );

    assert( config["fluid-solver"] );

    std::string fluidSolver = config["fluid-solver"].as<std::string>();

    assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" || fluidSolver == "compressible-solver" || fluidSolver == "sdc-pimple-solver" || fluidSolver == "sdc-laplacian-solver" || fluidSolver == "esdirk-pimple-solver" );

    std::shared_ptr<foamFluidSolver> fluid;
    std::shared_ptr<sdc::SDC> sdc;
    std::shared_ptr<sdc::ESDIRK> esdirk;

    if ( fluidSolver == "coupled-pressure-velocity-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "pimple-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "compressible-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( fluidSolver == "sdc-pimple-solver" || fluidSolver == "sdc-laplacian-solver" )
    {
        YAML::Node sdcConfig( config["sdc"] );
        assert( sdcConfig["convergence-tolerance"] );
        assert( sdcConfig["number-of-points"] );
        assert( sdcConfig["quadrature-rule"] );

        int n = sdcConfig["number-of-points"].as<int>();
        double tol = sdcConfig["convergence-tolerance"].as<double>();
        std::string quadratureRule = sdcConfig["quadrature-rule"].as<std::string>();

        std::shared_ptr<sdc::SDCSolver> solver;

        if ( fluidSolver == "sdc-pimple-solver" )
            solver = std::shared_ptr<sdc::SDCSolver>( new SDCFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

        if ( fluidSolver == "sdc-laplacian-solver" )
            solver = std::shared_ptr<sdc::SDCSolver>( new SDCLaplacianSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

        sdc = std::shared_ptr<sdc::SDC> ( new sdc::SDC( solver, quadratureRule, n, tol ) );
    }

    if ( fluidSolver == "esdirk-pimple-solver" )
    {
        YAML::Node esdirkConfig( config["esdirk"] );
        YAML::Node adaptiveTimeConfig( config["adaptive-time-stepping"] );

        assert( esdirkConfig["method"] );
        assert( adaptiveTimeConfig["enabled"] );

        std::string method = esdirkConfig["method"].as<std::string>();
        bool adaptiveTimeStepping = adaptiveTimeConfig["enabled"].as<bool>();
        std::string filter = "elementary";
        double adaptiveTolerance = 1.0e-3;
        double safetyFactor = 0.5;

        if ( adaptiveTimeStepping )
        {
            assert( adaptiveTimeConfig["filter"] );
            assert( adaptiveTimeConfig["tolerance"] );
            assert( adaptiveTimeConfig["safety-factor"] );
            filter = adaptiveTimeConfig["filter"].as<std::string>();
            adaptiveTolerance = adaptiveTimeConfig["tolerance"].as<double>();
            safetyFactor = adaptiveTimeConfig["safety-factor"].as<double>();
        }

        std::shared_ptr<sdc::SDCSolver> solver;
        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper;

        solver = std::shared_ptr<sdc::SDCSolver>( new SDCFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

        adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( adaptiveTimeStepping, filter, adaptiveTolerance, safetyFactor ) );

        esdirk = std::shared_ptr<sdc::ESDIRK>( new sdc::ESDIRK( solver, method, adaptiveTimeStepper ) );
    }

    assert( fluid || sdc || esdirk );

    if ( fluid )
        fluid->run();

    if ( sdc )
        sdc->run();

    if ( esdirk )
        esdirk->run();

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
