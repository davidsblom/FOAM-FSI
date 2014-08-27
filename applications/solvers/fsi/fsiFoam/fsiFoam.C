
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
#include "CoupledFluidSolver.H"
#include "FsiSolver.H"
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
#include "SpaceMappingSolver.H"

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

  if ( config["multi-level-acceleration"] )
  {
    std::string algorithm = config["multi-level-acceleration"]["algorithm"].as<std::string>();

    assert( !config["coupling-scheme-implicit"] );
    assert( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "ML-IQN-ILS" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" );

    int nbLevels = config["multi-level-acceleration"]["levels"].size();

    std::shared_ptr< std::deque<std::shared_ptr<SpaceMappingSolver> > > solvers;
    std::shared_ptr< std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > > models;

    solvers = std::shared_ptr< std::deque<std::shared_ptr<SpaceMappingSolver> > > ( new std::deque<std::shared_ptr<SpaceMappingSolver> >() );
    models = std::shared_ptr< std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> >  > ( new std::deque<std::shared_ptr<ImplicitMultiLevelFsiSolver> > () );

    // Create shared pointers to solvers
    std::shared_ptr<CoupledFluidSolver> fluid;
    std::shared_ptr<SolidSolver> solid;
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
    assert( configLevel["level"].as<int>() == level );
    assert( configLevel["extrapolation-order"] );
    assert( configLevel["max-iterations"] );
    assert( configPostProcessing["initial-relaxation"] );
    assert( configPostProcessing["timesteps-reused"] );
    assert( configPostProcessing["singularity-limit"] );
    assert( configPostProcessing["reuse-information-starting-from-time-index"] );

    bool parallel = config["parallel-coupling"].as<bool>();
    bool scaling = false;

    if ( parallel )
      scaling = true;

    int extrapolation = configLevel["extrapolation-order"].as<int>();
    int maxIter = configLevel["max-iterations"].as<int>();
    double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
    int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
    double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
    int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();

    assert( extrapolation >= 0 );
    assert( extrapolation <= 2 );
    assert( maxIter > 0 );
    assert( initialRelaxation > 0 );
    assert( initialRelaxation <= 1 );
    assert( nbReuse >= 0 );
    assert( singularityLimit < 1 );
    assert( singularityLimit > 0 );

    fluid = std::shared_ptr<CoupledFluidSolver> ( new CoupledFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );
    solid = std::shared_ptr<SolidSolver> ( new SolidSolver( "solid-level-" + std::to_string( level ), args, runTime ) );

    // Convergence measures
    convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

    setConvergenceMeasures( configMeasures, convergenceMeasures );

    multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, 0, nbLevels - 1 ) );
    multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, 1, nbLevels - 1 ) );

    int maxUsedIterations = solid->data.rows() * solid->data.cols();

    if ( parallel )
      maxUsedIterations += fluid->data.rows() * fluid->data.cols();

    multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    if ( algorithm == "output-space-mapping" )
      postProcessing = std::shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );
    else
      postProcessing = std::shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

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

      bool parallel = config["parallel-coupling"].as<bool>();
      int extrapolation = configLevel["extrapolation-order"].as<int>();
      int maxIter = configLevel["max-iterations"].as<int>();
      double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
      int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
      double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
      int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
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

      fluid = std::shared_ptr<CoupledFluidSolver> ( new CoupledFluidSolver( "fluid-level-" + std::to_string( level ), args, runTime ) );
      solid = std::shared_ptr<SolidSolver> ( new SolidSolver( "solid-level-" + std::to_string( level ), args, runTime ) );

      // Convergence measures
      convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

      setConvergenceMeasures( configMeasures, convergenceMeasures );

      multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModel->fsi->fluid, 0, level ) );
      multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModel->fsi->fluid, 1, level ) );

      multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

      maxUsedIterations = solid->data.rows() * solid->data.cols();

      if ( parallel )
        maxUsedIterations += fluid->data.rows() * fluid->data.cols();

      if ( algorithm == "output-space-mapping" )
        postProcessing = std::shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );
      else
        postProcessing = std::shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

      implicitMultiLevelFsiSolver = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

      models->push_back( implicitMultiLevelFsiSolver );
    }

    models->push_back( fineModel );

    for ( int level = 1; level < nbLevels; level++ )
    {
      YAML::Node configLevel( config["multi-level-acceleration"]["levels"][level] );
      YAML::Node configPostProcessing( config["multi-level-acceleration"]["levels"][level]["post-processing"] );

      if ( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" )
      {
        std::shared_ptr<SpaceMappingSolver> spaceMappingSolver;
        std::shared_ptr<SpaceMapping> spaceMapping;

        if ( level == 1 )
        {
          reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();

          if ( algorithm == "aggressive-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new AggressiveSpaceMapping( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

          if ( algorithm == "ASM-ILS" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ASMILS( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "manifold-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ManifoldMapping( models->at( level ), models->at( 0 ), maxIter, singularityLimit, nbReuse, reuseInformationStartingFromTimeIndex ) );

          if ( algorithm == "output-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new OutputSpaceMapping( models->at( level ), models->at( 0 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

          spaceMappingSolver = std::shared_ptr<SpaceMappingSolver>( new SpaceMappingSolver( models->at( level ), models->at( 0 ), spaceMapping ) );

          solvers->push_back( spaceMappingSolver );
        }

        if ( level > 1 )
        {
          reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();

          if ( algorithm == "aggressive-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new AggressiveSpaceMapping( models->at( level ), solvers->at( level - 1 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

          if ( algorithm == "ASM-ILS" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ASMILS( models->at( level ), solvers->at( level - 1 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

          if ( algorithm == "manifold-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new ManifoldMapping( models->at( level ), solvers->at( level - 2 ), maxIter, singularityLimit, nbReuse, reuseInformationStartingFromTimeIndex ) );

          if ( algorithm == "output-space-mapping" )
            spaceMapping = std::shared_ptr<SpaceMapping>( new OutputSpaceMapping( models->at( level ), solvers->at( level - 2 ), maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

          spaceMappingSolver = std::shared_ptr<SpaceMappingSolver>( new SpaceMappingSolver( models->at( level ), solvers->at( level - 2 ), spaceMapping ) );

          solvers->push_back( spaceMappingSolver );
        }
      }
    }

    std::shared_ptr<Solver> solver;

    if ( algorithm == "manifold-mapping" || algorithm == "output-space-mapping" || algorithm == "aggressive-space-mapping" || algorithm == "ASM-ILS" )
      solver = std::shared_ptr<Solver>( new MultiLevelSpaceMappingSolver( solvers, models ) );

    if ( algorithm == "ML-IQN-ILS" )
      solver = std::shared_ptr<Solver>( new MLIQNILSSolver( models ) );

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

    bool parallel = config["parallel-coupling"].as<bool>();
    int extrapolation = config["coupling-scheme-implicit"]["extrapolation-order"].as<int>();
    int maxIter = config["coupling-scheme-implicit"]["max-iterations"].as<int>();
    double initialRelaxation = configPostProcessing["initial-relaxation"].as<double>();
    int nbReuse = configPostProcessing["timesteps-reused"].as<int>();
    double singularityLimit = configPostProcessing["singularity-limit"].as<double>();
    int reuseInformationStartingFromTimeIndex = configPostProcessing["reuse-information-starting-from-time-index"].as<int>();
    bool scaling = false;
    std::string algorithm = configPostProcessing["algorithm"].as<std::string>();

    assert( algorithm == "IQN-ILS" or algorithm == "QN" or algorithm == "Aitken" );

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
    std::shared_ptr<CoupledFluidSolver> fluid;
    std::shared_ptr<SolidSolver> solid;
    std::shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
    std::shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
    std::shared_ptr<FsiSolver> fsi;
    std::shared_ptr< list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
    std::shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
    std::shared_ptr<PostProcessing> postProcessing;
    std::shared_ptr<ImplicitMultiLevelFsiSolver> implicitMultiLevelFsiSolver;

    fluid = std::shared_ptr<CoupledFluidSolver> ( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );
    solid = std::shared_ptr<SolidSolver> ( new SolidSolver( "solid", args, runTime ) );

    // Convergence measures
    convergenceMeasures = std::shared_ptr<list<std::shared_ptr<ConvergenceMeasure> > >( new list<std::shared_ptr<ConvergenceMeasure> > );

    setConvergenceMeasures( configMeasures, convergenceMeasures );

    multiLevelFluidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fluid, 0, 0 ) );
    multiLevelSolidSolver = std::shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fluid, 1, 0 ) );

    multiLevelFsiSolver = std::shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

    int maxUsedIterations = solid->data.rows() * solid->data.cols();

    if ( parallel )
      maxUsedIterations += fluid->data.rows() * fluid->data.cols();

    if ( algorithm == "Aitken" )
      postProcessing = std::shared_ptr<PostProcessing> ( new AitkenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter,  maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    if ( algorithm == "IQN-ILS" )
      postProcessing = std::shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

    if ( algorithm == "QN" )
      postProcessing = std::shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter,  maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

    implicitMultiLevelFsiSolver = std::shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

    implicitMultiLevelFsiSolver->run();
  }

  Info << "End\n" << endl;

  return (0);
}

// ************************************************************************* //
