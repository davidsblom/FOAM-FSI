
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <iostream>
#include <list>
#include <memory>
#include <yaml-cpp/yaml.h>

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
#include "SpaceMappingSolver.H"
#include "CompressibleFluidSolver.H"
#include "TPS.H"
#include "WendlandC0Function.H"
#include "WendlandC2Function.H"
#include "WendlandC4Function.H"
#include "WendlandC6Function.H"

using std::list;

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

            double tol = measure["relative-convergence-measure"]["limit"].as<double>();
            int dataId = measure["relative-convergence-measure"]["data-id"].as<int>();
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( dataId, tol ) ) );
        }

        if ( measure["min-iteration-convergence-measure"] )
        {
            assert( measure["min-iteration-convergence-measure"]["min-iterations"] );

            int minIter = measure["min-iteration-convergence-measure"]["min-iterations"].as<int>();
            convergenceMeasures->push_back( std::shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
        }
    }

    assert( convergenceMeasures->size() == configMeasures.size() );
    assert( convergenceMeasures->size() > 0 );
}

std::shared_ptr<rbf::RBFInterpolation> createRBFInterpolator(
    std::string interpolationFunction,
    double radius,
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

    std::shared_ptr<Time> runTime( new Time
        (
            Time::controlDictName,
            args->rootPath(),
            args->caseName()
        ) );

    // Load computation settings

    string filename = static_cast<std::string>( args->rootPath() ) + "/" + static_cast<std::string>( args->globalCaseName() ) + "/constant/fsi.yaml";

    YAML::Node config = YAML::LoadFile( filename );
    YAML::Node configInterpolation( config["interpolation"] );

    assert( config["fluid-solver"] );
    assert( config["solid-solver"] );

    std::string fluidSolver = config["fluid-solver"].as<std::string>();
    std::string solidSolver = config["solid-solver"].as<std::string>();

    assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" || fluidSolver == "compressible-solver" );
    assert( solidSolver == "segregated-solver" );

    assert( configInterpolation["coarsening"] );
    assert( configInterpolation["coarsening"]["enabled"] );
    assert( configInterpolation["radial-basis-function"] );
    assert( configInterpolation["radial-basis-function"]["function"] );
    bool coarsening = configInterpolation["coarsening"]["enabled"].as<bool>();
    std::string interpolationFunction = configInterpolation["radial-basis-function"]["function"].as<std::string>();
    double coarseningTol = 1.0e-5;
    int coarseningMinPoints = 1;
    int coarseningMaxPoints = 10;
    double radius = 1;
    bool livePointSelection = false;
    double tolLivePointSelection = 1.0e-5;
    bool cpu = false;

    assert( interpolationFunction == "thin-plate-spline" || interpolationFunction == "wendland-c0" || interpolationFunction == "wendland-c2" || interpolationFunction == "wendland-c4" || interpolationFunction == "wendland-c6" );

    if ( configInterpolation["radial-basis-function"]["cpu"] )
    {
        cpu = configInterpolation["radial-basis-function"]["cpu"].as<bool>();
    }

    if ( interpolationFunction != "thin-plate-spline" )
    {
        assert( configInterpolation["radial-basis-function"]["radius"] );
        radius = configInterpolation["radial-basis-function"]["radius"].as<double>();
    }

    if ( coarsening )
    {
        assert( configInterpolation["coarsening"]["tol"] );
        assert( configInterpolation["coarsening"]["min-points"] );
        assert( configInterpolation["coarsening"]["max-points"] );
        assert( configInterpolation["coarsening"]["live-point-selection"] );

        coarseningTol = configInterpolation["coarsening"]["tol"].as<double>();
        coarseningMinPoints = configInterpolation["coarsening"]["min-points"].as<int>();
        coarseningMaxPoints = configInterpolation["coarsening"]["max-points"].as<int>();
        livePointSelection = configInterpolation["coarsening"]["live-point-selection"].as<bool>();

        if ( livePointSelection )
        {
            assert( configInterpolation["coarsening"]["tol-live-point-selection"] );
            tolLivePointSelection = configInterpolation["coarsening"]["tol-live-point-selection"].as<double>();
        }
    }

    if ( config["multi-level-acceleration"] )
    {
        std::string algorithm = config["multi-level-acceleration"]["algorithm"].as<std::string>();

        assert( !config["coupling-scheme-implicit"] );
        assert( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "ML-IQN-ILS" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" );

        int nbLevels = config["multi-level-acceleration"]["levels"].size();

        std::shared_ptr< std::deque<std::shared_ptr<SpaceMappingSolver> > > solvers;
        std::shared_ptr< std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

        solvers = std::shared_ptr< std::deque<std::shared_ptr<SpaceMappingSolver> > > ( new std::deque<std::shared_ptr<SpaceMappingSolver> >() );
        models = std::shared_ptr< std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > > ( new std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > () );

        // Create shared pointers to solvers
        std::shared_ptr<foamFluidSolver> fluid;
        std::shared_ptr<foamSolidSolver> solid;
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
        double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
        int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
        double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
        int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
        int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
        double beta = configPostProcessing["beta"].as<double>();
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

        if ( solidSolver == "segregated-solver" )
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
            double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
            int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
            double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
            int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
            int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
            double beta = configPostProcessing["beta"].as<double>();
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

            if ( solidSolver == "segregated-solver" )
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

                int maxIter = configLevel["max-iterations"].as<int>();
                int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
                double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
                int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
                bool updateJacobian = configPostProcessing["update-jacobian"].as<bool>();
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
                    spaceMapping = shared_ptr<SpaceMapping> ( new ASMILS( fineModel, coarseModel, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

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
        assert( configPostProcessing["initial-relaxation"] );
        assert( configPostProcessing["timesteps-reused"] );
        assert( configPostProcessing["singularity-limit"] );
        assert( configPostProcessing["reuse-information-starting-from-time-index"] );
        assert( configPostProcessing["algorithm"] );
        assert( configPostProcessing["max-used-iterations"] );

        bool parallel = config["parallel-coupling"].as<bool>();
        int extrapolation = config["coupling-scheme-implicit"]["extrapolation-order"].as<int>();
        int maxIter = config["coupling-scheme-implicit"]["max-iterations"].as<int>();
        double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
        int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
        double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
        int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
        int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
        bool scaling = false;
        std::string algorithm = configPostProcessing["algorithm"].as<std::string>();
        double beta = 1;
        bool updateJacobian = false;

        assert( algorithm == "QN" or algorithm == "Aitken" or algorithm == "Anderson" );

        if ( algorithm == "Anderson" )
        {
            assert( configPostProcessing["beta"] );
            assert( configPostProcessing["update-jacobian"] );
            beta = configPostProcessing["beta"].as<double>();
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
        std::shared_ptr<foamSolidSolver> solid;
        std::shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
        std::shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
        std::shared_ptr<FsiSolver> fsi;
        std::shared_ptr< list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        std::shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
        std::shared_ptr<PostProcessing> postProcessing;
        std::shared_ptr<ImplicitMultiLevelFsiSolver> implicitMultiLevelFsiSolver;

        if ( fluidSolver == "coupled-pressure-velocity-solver" )
            fluid = std::shared_ptr<foamFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

        if ( fluidSolver == "pimple-solver" )
            fluid = std::shared_ptr<foamFluidSolver> ( new FluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

        if ( fluidSolver == "compressible-solver" )
            fluid = std::shared_ptr<foamFluidSolver> ( new CompressibleFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

        if ( solidSolver == "segregated-solver" )
        {
            std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );

            std::shared_ptr<rbf::RBFCoarsening> interpolator( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

            solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid", args, runTime, interpolator ) );
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

        multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 0, 0 ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToCouplingMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        rbfInterpolator = createRBFInterpolator( interpolationFunction, radius, cpu );
        rbfInterpToMesh = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, false, coarseningTol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, false ) );

        multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, rbfInterpToCouplingMesh, rbfInterpToMesh, 1, 0 ) );

        multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

        if ( algorithm == "Aitken" )
            postProcessing = std::shared_ptr<PostProcessing> ( new AitkenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

        if ( algorithm == "Anderson" )
            postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian ) );

        if ( algorithm == "QN" )
            postProcessing = std::shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex ) );

        implicitMultiLevelFsiSolver = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

        implicitMultiLevelFsiSolver->run();
    }

    Info << "End\n" << endl;

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    return (0);
}

// ************************************************************************* //
