
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <mpi.h>
#include <iostream>
#include <list>
#include <memory>
#include <yaml-cpp/yaml.h>

#include "dealiiSolidSolver.H"
#include "version.H"
#include "ASMILS.H"
#include "AggressiveSpaceMapping.H"
#include "AitkenPostProcessing.H"
#include "BroydenPostProcessing.H"
#include "ConvergenceMeasure.H"
#include "FluidSolver.H"
#include "CoupledFluidSolver.H"
#include "FsiSolver.H"
#include "AndersonPostProcessing.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MLIQNILSSolver.H"
#include "ManifoldMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MultiLevelFsiSolver.H"
#include "MultiLevelSolver.H"
#include "MultiLevelSpaceMappingSolver.H"
#include "OutputSpaceMapping.H"
#include "RBFCoarsening.H"
#include "RBFInterpolation.H"
#include "RelativeConvergenceMeasure.H"
#include "SolidSolver.H"
#include "ElasticSolidSolver.H"
#include "SpaceMappingSolver.H"
#include "CompressibleFluidSolver.H"
#include "TPSFunction.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"
#include "AdaptiveTimeStepper.H"
#include "SDCFsiSolverInterface.H"
#include "SDCDynamicMeshFluidSolver.H"
#include "SDCSolidSolver.H"
#include "SDCFsiSolver.H"
#include "PIES.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "AbsoluteConvergenceMeasure.H"
#include "SteadyStateFluidSolver.H"
#include "SteadyStateSolidSolver.H"

using std::list;

void setConvergenceMeasures(
    YAML::Node & configMeasures,
    std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > > & convergenceMeasures
    );

std::shared_ptr<rbf::RBFInterpolation> createRBFInterpolator(
    std::string interpolationFunction,
    scalar radius,
    bool cpu
    );

void setConvergenceMeasures(
    YAML::Node & configMeasures,
    std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > > & convergenceMeasures
    )
{
    for ( unsigned i = 0; i < configMeasures.size(); i++ )
    {
        YAML::Node measure( configMeasures[i] );

        if ( measure["relative-convergence-measure"] )
        {
            assert( measure["relative-convergence-measure"]["limit"] );
            assert( measure["relative-convergence-measure"]["data-id"] );
            assert( measure["relative-convergence-measure"]["suffices"] );

            scalar tol = measure["relative-convergence-measure"]["limit"].as<scalar>();
            int dataId = measure["relative-convergence-measure"]["data-id"].as<int>();
            bool suffices = measure["relative-convergence-measure"]["suffices"].as<bool>();
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( dataId, suffices, tol ) ) );
        }

        if ( measure["min-iteration-convergence-measure"] )
        {
            assert( measure["min-iteration-convergence-measure"]["min-iterations"] );
            assert( measure["min-iteration-convergence-measure"]["suffices"] );

            int minIter = measure["min-iteration-convergence-measure"]["min-iterations"].as<int>();
            bool suffices = measure["min-iteration-convergence-measure"]["suffices"].as<bool>();
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, suffices, minIter ) ) );
        }

        if ( measure["absolute-convergence-measure"] )
        {
            assert( measure["absolute-convergence-measure"]["limit"] );
            assert( measure["absolute-convergence-measure"]["data-id"] );
            assert( measure["absolute-convergence-measure"]["suffices"] );

            scalar tol = measure["absolute-convergence-measure"]["limit"].as<scalar>();
            int dataId = measure["absolute-convergence-measure"]["data-id"].as<int>();
            bool suffices = measure["absolute-convergence-measure"]["suffices"].as<bool>();
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new AbsoluteConvergenceMeasure( dataId, suffices, tol ) ) );
        }

        if ( measure["residual-relative-convergence-measure"] )
        {
            assert( measure["residual-relative-convergence-measure"]["limit"] );
            assert( measure["residual-relative-convergence-measure"]["data-id"] );
            assert( measure["residual-relative-convergence-measure"]["suffices"] );

            scalar tol = measure["residual-relative-convergence-measure"]["limit"].as<scalar>();
            int dataId = measure["residual-relative-convergence-measure"]["data-id"].as<int>();
            bool suffices = measure["residual-relative-convergence-measure"]["suffices"].as<bool>();
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new ResidualRelativeConvergenceMeasure( dataId, suffices, tol ) ) );
        }
    }

    assert( convergenceMeasures->size() == configMeasures.size() );
    assert( convergenceMeasures->size() > 0 );
}

std::shared_ptr<rbf::RBFInterpolation> createRBFInterpolator(
    std::string interpolationFunction,
    scalar radius,
    bool cpu
    )
{
    std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;

    if ( interpolationFunction == "thin-plate-spline" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::TPSFunction() );

    if ( interpolationFunction == "wendland-c0" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC0Function( radius ) );

    if ( interpolationFunction == "wendland-c2" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC2Function( radius ) );

    if ( interpolationFunction == "wendland-c4" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC4Function( radius ) );

    if ( interpolationFunction == "wendland-c6" )
        rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface>( new rbf::WendlandC6Function( radius ) );

    assert( rbfFunction );

    bool polynomialTerm = true;
    return std::shared_ptr<rbf::RBFInterpolation>( new rbf::RBFInterpolation( rbfFunction, polynomialTerm, cpu ) );
}

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

    Info << nl << "FOAM-FSI build: " << word( FOAMFSIbuild ) << nl << endl;

    std::shared_ptr<Time> runTime( new Time
        (
            Time::controlDictName,
            args->rootPath(),
            args->caseName()
        ) );

    // Load computation settings

    string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";

    YAML::Node config = YAML::LoadFile( filename );

    assert( config["interpolation"] );

    YAML::Node configInterpolation( config["interpolation"] );

    assert( config["fluid-solver"] );
    assert( config["solid-solver"] );

    std::string fluidSolver = config["fluid-solver"].as<std::string>();
    std::string solidSolver = config["solid-solver"].as<std::string>();

    assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" || fluidSolver == "compressible-solver" || fluidSolver == "steady-state-pimple-solver" );
    assert( solidSolver == "nonlinear-elastic-solver" || solidSolver == "dealii-solver" || solidSolver == "linear-elastic-solver" || solidSolver == "steady-state-nonlinear-elastic-solver" );

    assert( configInterpolation["coarsening"] );
    assert( configInterpolation["coarsening"]["enabled"] );
    assert( configInterpolation["radial-basis-function"] );
    assert( configInterpolation["radial-basis-function"]["function"] );
    bool coarsening = configInterpolation["coarsening"]["enabled"].as<bool>();
    std::string interpolationFunction = configInterpolation["radial-basis-function"]["function"].as<std::string>();
    scalar coarseningTol = 1.0e-5;
    int coarseningMinPoints = 1;
    int coarseningMaxPoints = 10;
    scalar radius = 1;
    bool livePointSelection = false;
    scalar tolLivePointSelection = 1.0e-5;
    bool cpu = false;

    assert( interpolationFunction == "thin-plate-spline" || interpolationFunction == "wendland-c0" || interpolationFunction == "wendland-c2" || interpolationFunction == "wendland-c4" || interpolationFunction == "wendland-c6" );

    if ( configInterpolation["radial-basis-function"]["cpu"] )
    {
        cpu = configInterpolation["radial-basis-function"]["cpu"].as<bool>();
    }

    if ( interpolationFunction != "thin-plate-spline" )
    {
        assert( configInterpolation["radial-basis-function"]["radius"] );
        radius = configInterpolation["radial-basis-function"]["radius"].as<scalar>();
    }

    if ( coarsening )
    {
        assert( configInterpolation["coarsening"]["tol"] );
        assert( configInterpolation["coarsening"]["min-points"] );
        assert( configInterpolation["coarsening"]["max-points"] );
        assert( configInterpolation["coarsening"]["live-point-selection"] );

        coarseningTol = configInterpolation["coarsening"]["tol"].as<scalar>();
        coarseningMinPoints = configInterpolation["coarsening"]["min-points"].as<int>();
        coarseningMaxPoints = configInterpolation["coarsening"]["max-points"].as<int>();
        livePointSelection = configInterpolation["coarsening"]["live-point-selection"].as<bool>();

        if ( livePointSelection )
        {
            assert( configInterpolation["coarsening"]["tol-live-point-selection"] );
            tolLivePointSelection = configInterpolation["coarsening"]["tol-live-point-selection"].as<scalar>();
        }
    }

    if ( config["multi-level-acceleration"] )
    {
        std::string algorithm = config["multi-level-acceleration"]["algorithm"].as<std::string>();

        assert( !config["coupling-scheme-implicit"] );
        assert( solidSolver == "nonlinear-elastic-solver" );
        assert( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "ML-IQN-ILS" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" );

        int nbLevels = config["multi-level-acceleration"]["levels"].size();

        std::shared_ptr< std::deque<std::shared_ptr<SpaceMappingSolver> > > solvers;
        std::shared_ptr< std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

        solvers = std::shared_ptr< std::deque<std::shared_ptr<SpaceMappingSolver> > > ( new std::deque<std::shared_ptr<SpaceMappingSolver> >() );
        models = std::shared_ptr< std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > > ( new std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > () );

        // Create shared pointers to solvers
        std::shared_ptr<foamFluidSolver> fluid;
        std::shared_ptr<BaseMultiLevelSolver> solid;
        std::shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        std::shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        std::shared_ptr<FsiSolver> fsi;
        std::shared_ptr< list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        std::shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        std::shared_ptr<PostProcessing> postProcessing;
        std::shared_ptr<ImplicitMultiLevelFsiSolver> implicitMultiLevelFsiSolver;
        std::shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

        int level = nbLevels - 1;

        YAML::Node configLevel( config["multi-level-acceleration"]["levels"][level] );
        YAML::Node configPostProcessing( config["multi-level-acceleration"]["levels"][level]["post-processing"] );
        YAML::Node configMeasures( config["multi-level-acceleration"]["levels"][level]["convergence-measures"] );

        assert( config["parallel-coupling"] );
        assert( config["multi-level-acceleration"]["synchronization"] );
        assert( configLevel["level"].as<int>() == level );
        assert( configLevel["extrapolation-order"] );
        assert( configLevel["max-iterations"] );
        assert( configPostProcessing["initial-relaxation"] );
        assert( configPostProcessing["timesteps-reused"] );
        assert( configPostProcessing["singularity-limit"] );
        assert( configPostProcessing["reuse-information-starting-from-time-index"] );
        assert( configPostProcessing["max-used-iterations"] );
        assert( configPostProcessing["beta"] );
        assert( configPostProcessing["update-jacobian"] );

        int order = 0;

        if ( algorithm == "output-space-mapping" )
        {
            assert( config["multi-level-acceleration"]["order"] );
            order = config["multi-level-acceleration"]["order"].as<int>();
        }

        bool parallel = config["parallel-coupling"].as<bool>();
        bool scaling = false;
        bool synchronization = config["multi-level-acceleration"]["synchronization"].as<bool>();

        if ( parallel )
            scaling = true;

        int extrapolation = configLevel["extrapolation-order"].as<int>();
        int maxIter = configLevel["max-iterations"].as<int>();
        scalar initialRelaxation = configPostProcessing["initial-relaxation"].as<scalar>();
        int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
        scalar singularityLimit = configPostProcessing["singularity-limit"].as<scalar>();
        int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
        int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
        scalar beta = configPostProcessing["beta"].as<scalar>();
        bool updateJacobian = configPostProcessing["update-jacobian"].as<bool>();

        assert( extrapolation >= 0 );
        assert( extrapolation <= 2 );
        assert( maxIter > 0 );
        assert( initialRelaxation > 0 );
        assert( initialRelaxation <= 1 );
        assert( nbReuse >= 0 );
        assert( singularityLimit < 1 );
        assert( singularityLimit > 0 );

        if ( fluidSolver == "coupled-pressure-velocity-solver" )
            fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

        if ( fluidSolver == "pimple-solver" )
            fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

        if ( fluidSolver == "compressible-solver" )
            fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

        if ( solidSolver == "nonlinear-elastic-solver" )
        {
            std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );

            std::shared_ptr<rbf::RBFCoarsening> interpolator( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

            solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid-level-" + std::to_string( level ), args, runTime, interpolator ) );
        }

        // Convergence measures
        convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

        setConvergenceMeasures( configMeasures, convergenceMeasures );

        std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
        std::shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
        std::shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, nbLevels - 1 ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, nbLevels - 1 ) );

        multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        fineModel = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        for ( int level = 0; level < nbLevels - 1; level++ )
        {
            YAML::Node configLevel( config["multi-level-acceleration"]["levels"][level] );
            YAML::Node configPostProcessing( config["multi-level-acceleration"]["levels"][level]["post-processing"] );
            YAML::Node configMeasures( config["multi-level-acceleration"]["levels"][level]["convergence-measures"] );

            assert( configLevel["level"].as<int>() == level );
            assert( configLevel["extrapolation-order"] );
            assert( configLevel["max-iterations"] );
            assert( configPostProcessing["initial-relaxation"] );
            assert( configPostProcessing["timesteps-reused"] );
            assert( configPostProcessing["singularity-limit"] );
            assert( configPostProcessing["reuse-information-starting-from-time-index"] );
            assert( configPostProcessing["max-used-iterations"] );
            assert( configPostProcessing["beta"] );
            assert( configPostProcessing["update-jacobian"] );

            bool parallel = config["parallel-coupling"].as<bool>();
            int extrapolation = configLevel["extrapolation-order"].as<int>();
            int maxIter = configLevel["max-iterations"].as<int>();
            scalar initialRelaxation = configPostProcessing["initial-relaxation"].as<scalar>();
            int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
            scalar singularityLimit = configPostProcessing["singularity-limit"].as<scalar>();
            int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
            int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
            scalar beta = configPostProcessing["beta"].as<scalar>();
            bool scaling = false;
            bool updateJacobian = configPostProcessing["update-jacobian"].as<bool>();

            if ( parallel )
                scaling = true;

            assert( extrapolation >= 0 );
            assert( extrapolation <= 2 );
            assert( maxIter > 0 );
            assert( initialRelaxation > 0 );
            assert( initialRelaxation <= 1 );
            assert( nbReuse >= 0 );
            assert( singularityLimit < 1 );
            assert( singularityLimit > 0 );

            if ( fluidSolver == "coupled-pressure-velocity-solver" )
                fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

            if ( fluidSolver == "pimple-solver" )
                fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

            if ( fluidSolver == "compressible-solver" )
                fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

            if ( solidSolver == "nonlinear-elastic-solver" )
            {
                std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );

                std::shared_ptr<rbf::RBFCoarsening> interpolator( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

                solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid-level-" + std::to_string( level ), args, runTime, interpolator ) );
            }

            // Convergence measures
            convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

            setConvergenceMeasures( configMeasures, convergenceMeasures );

            std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
            std::shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
            std::shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

            rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
            rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

            rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
            rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

            multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, level ) );

            rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
            rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

            rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
            rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

            multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, level ) );

            multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

            postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

            implicitMultiLevelFsiSolver = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

            models->push_back( implicitMultiLevelFsiSolver );
        }

        models->push_back( fineModel );

        if ( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" )
        {
            for ( int level = 0; level < nbLevels - 1; level++ )
            {
                YAML::Node configLevel( config["multi-level-acceleration"]["levels"][level] );
                YAML::Node configPostProcessing( configLevel["post-processing"] );

                assert( config["multi-level-acceleration"]["initial-solution-coarse-model"] );
                assert( configPostProcessing["beta"] );

                int maxIter = configLevel["max-iterations"].as<int>();
                int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
                scalar singularityLimit = configPostProcessing["singularity-limit"].as<scalar>();
                int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
                bool updateJacobian = configPostProcessing["update-jacobian"].as<bool>();
                scalar beta = configPostProcessing["beta"].as<scalar>();
                bool initialSolutionCoarseModel = config["multi-level-acceleration"]["initial-solution-coarse-model"].as<bool>();

                shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
                shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel;

                fineModel = models->at( level + 1 );

                if ( level == 0 )
                    coarseModel = models->at( level );

                if ( level > 0 )
                    coarseModel = solvers->at( level - 1 );

                shared_ptr<SpaceMapping> spaceMapping;

                if ( algorithm == "manifold-mapping" )
                    spaceMapping = shared_ptr<SpaceMapping> ( new ManifoldMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian, initialSolutionCoarseModel ) );

                if ( algorithm == "output-space-mapping" )
                    spaceMapping = shared_ptr<SpaceMapping> ( new OutputSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

                if ( algorithm == "aggressive-space-mapping" )
                    spaceMapping = shared_ptr<SpaceMapping> ( new AggressiveSpaceMapping( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

                if ( algorithm == "ASM-ILS" )
                    spaceMapping = shared_ptr<SpaceMapping> ( new ASMILS( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, beta ) );

                shared_ptr<SpaceMappingSolver > spaceMappingSolver( new SpaceMappingSolver( fineModel, coarseModel, spaceMapping ) );

                solvers->push_back( spaceMappingSolver );

                spaceMapping.reset();
                spaceMappingSolver.reset();
            }
        }

        std::shared_ptr<Solver> solver;

        if ( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" )
            solver = std::shared_ptr<Solver>( new MultiLevelSpaceMappingSolver( solvers, models, synchronization ) );

        if ( algorithm == "ML-IQN-ILS" )
            solver = std::shared_ptr<Solver>( new MLIQNILSSolver( models, synchronization ) );

        solver->run();
    }
    else
    {
        YAML::Node configPostProcessing( config["coupling-scheme-implicit"]["post-processing"] );
        YAML::Node configMeasures( config["coupling-scheme-implicit"]["convergence-measures"] );

        assert( config["parallel-coupling"] );
        assert( config["coupling-scheme-implicit"]["extrapolation-order"] );
        assert( config["coupling-scheme-implicit"]["max-iterations"] );
        assert( config["time-integration-scheme"] );
        assert( config["adaptive-time-stepping"] );
        assert( configPostProcessing["initial-relaxation"] );
        assert( configPostProcessing["timesteps-reused"] );
        assert( configPostProcessing["singularity-limit"] );
        assert( configPostProcessing["reuse-information-starting-from-time-index"] );
        assert( configPostProcessing["algorithm"] );
        assert( configPostProcessing["max-used-iterations"] );

        bool parallel = config["parallel-coupling"].as<bool>();
        int extrapolation = config["coupling-scheme-implicit"]["extrapolation-order"].as<int>();
        int maxIter = config["coupling-scheme-implicit"]["max-iterations"].as<int>();
        scalar initialRelaxation = configPostProcessing["initial-relaxation"].as<scalar>();
        int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
        scalar singularityLimit = configPostProcessing["singularity-limit"].as<scalar>();
        int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
        int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
        bool scaling = false;
        std::string algorithm = configPostProcessing["algorithm"].as<std::string>();
        scalar beta = 1;
        bool updateJacobian = false;
        std::string firstParticipant = "fluid-solver";
        std::string timeIntegrationScheme = config["time-integration-scheme"].as<std::string>();

        if ( not parallel )
        {
            assert( config["first-participant"] );
            firstParticipant = config["first-participant"].as<std::string>();
            assert( firstParticipant == "fluid-solver" || firstParticipant == "solid-solver" );
        }

        assert( algorithm == "QN" or algorithm == "Aitken" or algorithm == "Anderson" );
        assert( timeIntegrationScheme == "bdf" || timeIntegrationScheme == "esdirk" || timeIntegrationScheme == "sdc" || timeIntegrationScheme == "picard-integral-exponential-solver" );

        std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper;

        YAML::Node adaptiveTimeConfig( config["adaptive-time-stepping"] );
        assert( adaptiveTimeConfig["enabled"] );
        bool adaptiveTimeStepping = adaptiveTimeConfig["enabled"].as<bool>();

        std::string filter = "elementary";
        scalar adaptiveTolerance = 1.0e-3;
        scalar safetyFactor = 0.5;

        if ( adaptiveTimeStepping )
        {
            assert( adaptiveTimeConfig["filter"] );
            assert( adaptiveTimeConfig["tolerance"] );
            assert( adaptiveTimeConfig["safety-factor"] );
            filter = adaptiveTimeConfig["filter"].as<std::string>();
            adaptiveTolerance = adaptiveTimeConfig["tolerance"].as<scalar>();
            safetyFactor = adaptiveTimeConfig["safety-factor"].as<scalar>();
        }

        adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> ( new sdc::AdaptiveTimeStepper( adaptiveTimeStepping, filter, adaptiveTolerance, safetyFactor ) );

        if ( algorithm == "Anderson" )
        {
            assert( configPostProcessing["beta"] );
            assert( configPostProcessing["update-jacobian"] );
            beta = configPostProcessing["beta"].as<scalar>();
            updateJacobian = configPostProcessing["update-jacobian"].as<bool>();
        }

        if ( parallel )
            scaling = true;

        assert( extrapolation >= 0 );
        assert( extrapolation <= 2 );
        assert( maxIter > 0 );
        assert( initialRelaxation > 0 );
        assert( initialRelaxation <= 1 );
        assert( nbReuse >= 0 );
        assert( singularityLimit < 1 );
        assert( singularityLimit > 0 );

        // Create shared pointers to solvers
        std::shared_ptr<foamFluidSolver> fluid;
        std::shared_ptr<BaseMultiLevelSolver> solid;
        std::shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        std::shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        std::shared_ptr<FsiSolver> fsi;
        std::shared_ptr< list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        std::shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        std::shared_ptr<PostProcessing> postProcessing;
        std::shared_ptr<ImplicitMultiLevelFsiSolver> implicitMultiLevelFsiSolver;

        std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver;
        std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver;

        if ( fluidSolver == "coupled-pressure-velocity-solver" )
        {
            assert( timeIntegrationScheme == "bdf" );
            assert( not adaptiveTimeStepping );

            fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
        }

        if ( fluidSolver == "pimple-solver" )
        {
            if ( timeIntegrationScheme == "bdf" )
            {
                assert( not adaptiveTimeStepping );

                fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
            }

            if ( timeIntegrationScheme == "esdirk" || timeIntegrationScheme == "sdc" || timeIntegrationScheme == "picard-integral-exponential-solver" )
                sdcFluidSolver = std::shared_ptr<sdc::SDCFsiSolverInterface> ( new SDCDynamicMeshFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
        }

        if ( fluidSolver == "steady-state-pimple-solver" )
        {
            assert( not adaptiveTimeStepping );

            fluid = std::shared_ptr<foamFluidSolver> ( new SteadyStateFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
        }

        if ( fluidSolver == "compressible-solver" )
        {
            assert( timeIntegrationScheme == "bdf" );
            assert( not adaptiveTimeStepping );

            fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
        }

        if ( solidSolver == "nonlinear-elastic-solver" )
        {
            if ( timeIntegrationScheme == "bdf" )
            {
                assert( not adaptiveTimeStepping );

                solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid", args, runTime ) );
            }

            if ( timeIntegrationScheme == "esdirk" || timeIntegrationScheme == "sdc" || timeIntegrationScheme == "picard-integral-exponential-solver" )
                sdcSolidSolver = std::shared_ptr<sdc::SDCFsiSolverInterface> ( new SDCSolidSolver( "solid", args, runTime ) );
        }

        if ( solidSolver == "steady-state-nonlinear-elastic-solver" )
        {
            assert( not adaptiveTimeStepping );

            solid = std::shared_ptr<foamSolidSolver> ( new SteadyStateSolidSolver( "solid", args, runTime ) );
        }

        if ( solidSolver == "linear-elastic-solver" )
        {
            assert( not adaptiveTimeStepping );
            assert( timeIntegrationScheme == "bdf" );

            solid = std::shared_ptr<foamSolidSolver> ( new ElasticSolidSolver( "solid", args, runTime ) );
        }

        if ( solidSolver == "dealii-solver" )
        {
            assert( not adaptiveTimeStepping );

            dealiifsi::DataStorage data;
            data.read_data( "deal-fsi.prm" );
            data.time_step = runTime->deltaT().value();
            data.final_time = runTime->endTime().value();

            if ( timeIntegrationScheme == "bdf" )
            {
                assert( not adaptiveTimeStepping );

                solid = std::shared_ptr<BaseMultiLevelSolver>( new dealiiSolidSolver<2> ( data ) );
            }

            if ( timeIntegrationScheme == "esdirk" || timeIntegrationScheme == "sdc" || timeIntegrationScheme == "picard-integral-exponential-solver" )
                sdcSolidSolver = std::shared_ptr<sdc::SDCFsiSolverInterface> ( new dealiiSolidSolver<2> ( data ) );
        }

        if ( timeIntegrationScheme == "esdirk" || timeIntegrationScheme == "sdc" || timeIntegrationScheme == "picard-integral-exponential-solver" )
        {
            assert( sdcFluidSolver );
            assert( sdcSolidSolver );
            assert( not fluid );
            assert( not solid );
            fluid = std::dynamic_pointer_cast<foamFluidSolver>( sdcFluidSolver );
            solid = std::dynamic_pointer_cast<BaseMultiLevelSolver>( sdcSolidSolver );
        }

        assert( fluid );
        assert( solid );

        // Convergence measures
        convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

        setConvergenceMeasures( configMeasures, convergenceMeasures );

        std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator;
        std::shared_ptr<rbf::RBFCoarsening> rbfInterpToCouplingMesh;
        std::shared_ptr<rbf::RBFCoarsening> rbfInterpToMesh;

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        if ( firstParticipant == "fluid-solver" )
            multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        if ( firstParticipant == "solid-solver" )
            multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 1 ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        if ( firstParticipant == "fluid-solver" )
            multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        if ( firstParticipant == "solid-solver" )
            multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 1 ) );

        if ( firstParticipant == "fluid-solver" )
            multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        if ( firstParticipant == "solid-solver" )
            multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelSolidSolver, multiLevelFluidSolver, convergenceMeasures, parallel, extrapolation ) );

        if ( algorithm == "Aitken" )
            postProcessing = std::shared_ptr<PostProcessing> ( new AitkenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

        if ( algorithm == "Anderson" )
            postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        if ( algorithm == "QN" )
            postProcessing = std::shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex ) );

        if ( timeIntegrationScheme == "bdf" )
        {
            implicitMultiLevelFsiSolver = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );
            implicitMultiLevelFsiSolver->run();
        }

        if ( timeIntegrationScheme == "esdirk" || timeIntegrationScheme == "sdc" || timeIntegrationScheme == "picard-integral-exponential-solver" )
        {
            std::shared_ptr<SDCFsiSolver> sdcFsiSolver( new SDCFsiSolver( sdcFluidSolver, sdcSolidSolver, postProcessing ) );
            std::shared_ptr<sdc::TimeIntegrationScheme> timeSolver;

            if ( timeIntegrationScheme == "esdirk" )
            {
                YAML::Node esdirkConfig( config["esdirk"] );

                assert( esdirkConfig["method"] );
                assert( adaptiveTimeStepper );

                std::string method = esdirkConfig["method"].as<std::string>();

                timeSolver = std::shared_ptr<sdc::TimeIntegrationScheme>( new sdc::ESDIRK( sdcFsiSolver, method, adaptiveTimeStepper ) );
            }

            if ( timeIntegrationScheme == "sdc" )
            {
                assert( config["sdc"] );
                YAML::Node sdcConfig( config["sdc"] );
                assert( sdcConfig["convergence-tolerance"] );
                assert( sdcConfig["number-of-points"] );
                assert( sdcConfig["quadrature-rule"] );
                assert( sdcConfig["min-sweeps"] );
                assert( sdcConfig["max-sweeps"] );
                assert( adaptiveTimeStepper );

                int n = sdcConfig["number-of-points"].as<int>();
                scalar tol = sdcConfig["convergence-tolerance"].as<scalar>();
                std::string quadratureRule = sdcConfig["quadrature-rule"].as<std::string>();
                int minSweeps = sdcConfig["min-sweeps"].as<int>();
                int maxSweeps = sdcConfig["max-sweeps"].as<int>();

                timeSolver = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( sdcFsiSolver, adaptiveTimeStepper, quadratureRule, n, tol, minSweeps, maxSweeps ) );
            }

            if ( timeIntegrationScheme == "picard-integral-exponential-solver" )
            {
                assert( config["picard-integral-exponential-solver"] );
                YAML::Node piesConfig( config["picard-integral-exponential-solver"] );
                assert( piesConfig["delta"] );
                assert( piesConfig["convergence-tolerance"] );
                assert( piesConfig["min-sweeps"] );
                assert( piesConfig["max-sweeps"] );
                assert( piesConfig["rho"] );

                scalar delta = piesConfig["delta"].as<scalar>();
                scalar tol = piesConfig["convergence-tolerance"].as<scalar>();
                scalar rho = piesConfig["rho"].as<scalar>();
                int minSweeps = piesConfig["min-sweeps"].as<int>();
                int maxSweeps = piesConfig["max-sweeps"].as<int>();

                timeSolver = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::PIES( sdcFsiSolver, adaptiveTimeStepper, rho, delta, tol, minSweeps, maxSweeps ) );
            }

            assert( timeSolver );

            timeSolver->run();
        }
    }

    Info << "End\n" << endl;

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    return (0);
}

// ************************************************************************* //
