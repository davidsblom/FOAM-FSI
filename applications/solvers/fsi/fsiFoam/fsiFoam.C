
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
#include "steadyCoupledFluidSolver.H"
#include "FsiSolver.H"
#include "AndersonPostProcessing.H"
#include "IQNILSPostProcessing.H"
#include "ImplicitMultiLevelFsiSolver.H"
#include "MLIQNILSSolver.H"
#include "ManifoldMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MultiLevelFsiSolver.H"
#include "MultiLevelSolver.H"
#include "MultiLevelSpaceMappingSolver.H"
#include "OutputSpaceMapping.H"
#include "RBFInterpolation.H"
#include "RelativeConvergenceMeasure.H"
#include "SolidSolver.H"
#include "steadySolidSolver.H"
#include "SpaceMappingSolver.H"
#include "ResponseParameterMapping.H"

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

  assert( config["fluid-solver"] );

  std::string fluidSolver = config["fluid-solver"].as<std::string>();
  std::string solidSolver = config["solid-solver"].as<std::string>();

  assert( fluidSolver == "coupled-pressure-velocity-solver" || fluidSolver == "pimple-solver" || fluidSolver == "steady-state-coupled-pressure-velocity-solver" );
  assert( solidSolver == "segregated-solver" || solidSolver == "steady-state-segregated-solver" );

  if ( config["multi-level-acceleration"] )
  {
    std::string algorithm = config["multi-level-acceleration"]["algorithm"].as<std::string>();

    assert( !config["coupling-scheme-implicit"] );
    assert( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "ML-IQN-ILS" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" || algorithm == "response-parameter-mapping" );

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
    double beta = configPostProcessing["beta"].as<int>();

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

    if ( fluidSolver == "steady-state-coupled-pressure-velocity-solver" )
      fluid = std::shared_ptr<foamFluidSolver> ( new steadyCoupledFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

    if ( solidSolver == "segregated-solver" )
      solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid-level-" + std::to_string( level ), args, runTime ) );

    if ( solidSolver == "steady-state-segregated-solver" )
      solid = std::shared_ptr<foamSolidSolver> ( new steadySolidSolver( "solid-level-" + std::to_string( level ), args, runTime ) );

    // Convergence measures
    convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

    setConvergenceMeasures( configMeasures, convergenceMeasures );

    multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, 0, nbLevels - 1 ) );
    multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, 1, nbLevels - 1 ) );

    multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta ) );

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

      bool parallel = config["parallel-coupling"].as<bool>();
      int extrapolation = configLevel["extrapolation-order"].as<int>();
      int maxIter = configLevel["max-iterations"].as<int>();
      double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
      int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
      double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
      int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
      int maxUsedIterations = configPostProcessing["max-used-iterations"].as<int>();
      double beta = configPostProcessing["beta"].as<int>();
      bool scaling = false;

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

      if ( fluidSolver == "steady-state-coupled-pressure-velocity-solver" )
        fluid = std::shared_ptr<foamFluidSolver> ( new steadyCoupledFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );

      if ( solidSolver == "segregated-solver" )
        solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid-level-" + std::to_string( level ), args, runTime ) );

      if ( solidSolver == "steady-state-segregated-solver" )
        solid = std::shared_ptr<foamSolidSolver> ( new steadySolidSolver( "solid-level-" + std::to_string( level ), args, runTime ) );

      // Convergence measures
      convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

      setConvergenceMeasures( configMeasures, convergenceMeasures );

      multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, 0, level ) );
      multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->fluid, 1, level ) );

      multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

      postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta ) );

      implicitMultiLevelFsiSolver = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

      models->push_back( implicitMultiLevelFsiSolver );
    }

    models->push_back( fineModel );

    for ( int level = 1; level < nbLevels; level++ )
    {
      YAML::Node configLevel( config["multi-level-acceleration"]["levels"][level] );
      YAML::Node configPostProcessing( config["multi-level-acceleration"]["levels"][level]["post-processing"] );

      if ( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" || algorithm == "response-parameter-mapping" )
      {
        std::shared_ptr<SpaceMappingSolver> spaceMappingSolver;
        std::shared_ptr<SpaceMapping> spaceMapping;

        if ( level == 1 )
        {
          reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();

          if ( algorithm == "aggressive-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new AggressiveSpaceMapping( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "ASM-ILS" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ASMILS( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "manifold-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ManifoldMapping( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "response-parameter-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ResponseParameterMapping( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "output-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new OutputSpaceMapping( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

          spaceMappingSolver = std::shared_ptr<SpaceMappingSolver>( new SpaceMappingSolver( models->at( level ), models->at( 0 ), spaceMapping ) );

          solvers->push_back( spaceMappingSolver );
        }

        if ( level > 1 )
        {
          reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();

          if ( algorithm == "aggressive-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new AggressiveSpaceMapping( models->at( level ), solvers->at( level - 1 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "ASM-ILS" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ASMILS( models->at( level ), solvers->at( level - 1 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "manifold-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ManifoldMapping( models->at( level ), solvers->at( level - 2 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "response-parameter-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ResponseParameterMapping( models->at( level ), solvers->at( level - 2 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "output-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new OutputSpaceMapping( models->at( level ), solvers->at( level - 2 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, order ) );

          spaceMappingSolver = std::shared_ptr<SpaceMappingSolver>( new SpaceMappingSolver( models->at( level ), solvers->at( level - 2 ), spaceMapping ) );

          solvers->push_back( spaceMappingSolver );
        }
      }
    }

    std::shared_ptr<Solver> solver;

    if ( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" || algorithm == "response-parameter-mapping" )
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

    assert( algorithm == "IQN-ILS" or algorithm == "QN" or algorithm == "Aitken" or algorithm == "Anderson" );

    if ( algorithm == "Anderson" )
    {
      assert( configPostProcessing["beta"] );
      beta = configPostProcessing["beta"].as<double>();
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

    if ( fluidSolver == "steady-state-coupled-pressure-velocity-solver" )
      fluid = std::shared_ptr<foamFluidSolver> ( new steadyCoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    if ( solidSolver == "segregated-solver" )
      solid = std::shared_ptr<foamSolidSolver> ( new SolidSolver( "solid", args, runTime ) );

    if ( solidSolver == "steady-state-segregated-solver" )
      solid = std::shared_ptr<foamSolidSolver> ( new steadySolidSolver( "solid", args, runTime ) );

    // Convergence measures
    convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

    setConvergenceMeasures( configMeasures, convergenceMeasures );

    multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, 0, 0 ) );
    multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, 1, 0 ) );

    multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    if ( algorithm == "Aitken" )
      postProcessing = std::shared_ptr<PostProcessing> ( new AitkenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    if ( algorithm == "IQN-ILS" )
      postProcessing = std::shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    if ( algorithm == "Anderson" )
      postProcessing = std::shared_ptr<PostProcessing> ( new AndersonPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta ) );

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
