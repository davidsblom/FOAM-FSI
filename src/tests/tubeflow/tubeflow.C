
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include "ASMILS.H"
#include "AggressiveSpaceMapping.H"
#include "BroydenPostProcessing.H"
#include "IQNILSPostProcessing.H"
#include "MLIQNILSSolver.H"
#include "ManifoldMapping.H"
#include "MinIterationConvergenceMeasure.H"
#include "MultiLevelSpaceMappingSolver.H"
#include "OutputSpaceMapping.H"
#include "RelativeConvergenceMeasure.H"
#include "TubeFlowFluidSolver.H"
#include "TubeFlowSolidSolver.H"

using namespace std;
using namespace fsi;
using namespace tubeflow;

int main(
  int argc,
  char ** argv
  )
{
  // List of solvers
  list<std::string> fsiSolvers;
  fsiSolvers.push_back( "IQN-ILS" );
  fsiSolvers.push_back( "IQN-Broyden" );
  fsiSolvers.push_back( "ML-IQN-ILS" );
  fsiSolvers.push_back( "MM" );
  fsiSolvers.push_back( "OSM" );
  fsiSolvers.push_back( "ASM" );

  // fsiSolvers.push_back( "ASM-ILS" );

  // Physical settings
  double r0 = 0.2;
  double a0 = M_PI * r0 * r0;
  double u0 = 0.1;
  double p0 = 0;
  double dt = 0.1;
  double L = 1;
  double T = 10;
  double rho = 1.225;
  double E = 490;
  double h = 1.0e-3;
  double cmk = std::sqrt( E * h / (2 * rho * r0) );

  // Computational settings
  bool parallel;
  double tol = 1.0e-5;
  int maxIter = 100;
  double initialRelaxation = 1.0e-3;
  double singularityLimit = 1.0e-11;
  int extrapolation = 0;
  int minIter = 1;
  int reuseInformationStartingFromTimeIndex = 0;
  bool scaling = false;

  // Grid settings
  int N = pow( 10, 3 );
  int couplingGridSize = N;

  for ( int iParallel = 0; iParallel < 2; iParallel++ )
  {
    for ( int iReuse = 0; iReuse < 3; iReuse++ )
    {
      int nbReuse = 0;

      if ( iReuse == 1 )
        nbReuse = 4;

      if ( iReuse == 2 )
        nbReuse = 8;

      for ( int nbLevels = 2; nbLevels < 3; nbLevels++ )
      {
        for ( list<std::string>::iterator fsiSolver = fsiSolvers.begin(); fsiSolver != fsiSolvers.end(); ++fsiSolver )
        {
          if ( iParallel == 0 )
            parallel = false;

          if ( iParallel == 1 )
            parallel = true;

          if ( (*fsiSolver == "IQN-ILS" || *fsiSolver == "IQN-Broyden") && nbLevels == 3 )
            continue;

          std::string label = *fsiSolver + "_" + to_string( nbReuse ) + "_" + to_string( parallel );
          label += "_" + to_string( nbLevels );

          ifstream ifile( label + ".log" );

          if ( ifile )
            continue;            // log file exists;

          // Create shared pointers to solvers
          shared_ptr<TubeFlowFluidSolver> fluid;
          shared_ptr<TubeFlowSolidSolver> solid;
          shared_ptr<TubeFlowFluidSolver> fineModelFluid;
          shared_ptr<TubeFlowSolidSolver> fineModelSolid;
          shared_ptr<MultiLevelSolver> multiLevelFluidSolver;
          shared_ptr<MultiLevelSolver> multiLevelSolidSolver;
          shared_ptr<MultiLevelFsiSolver> fsi;
          shared_ptr< list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
          shared_ptr<MultiLevelFsiSolver> multiLevelFsiSolver;
          shared_ptr<PostProcessing> postProcessing;
          shared_ptr< deque<shared_ptr<ImplicitMultiLevelFsiSolver> > > models;
          shared_ptr< deque<shared_ptr<SpaceMappingSolver> > > solvers;
          shared_ptr<Solver> solver;
          shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;

          models = shared_ptr< deque<shared_ptr<ImplicitMultiLevelFsiSolver> >  > ( new deque<shared_ptr<ImplicitMultiLevelFsiSolver> > () );
          solvers = shared_ptr< deque<shared_ptr<SpaceMappingSolver> > > ( new deque<shared_ptr<SpaceMappingSolver> >() );

          if ( *fsiSolver == "MM" || *fsiSolver == "OSM" || *fsiSolver == "ASM" || *fsiSolver == "ML-IQN-ILS" || *fsiSolver == "ASM-ILS" )
          {
            minIter = 1;
            tol = 1.0e-5;
            N = couplingGridSize;

            fineModelFluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
            fineModelSolid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
            multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelFluid, fineModelFluid, 0, nbLevels - 1 ) );
            multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fineModelSolid, fineModelSolid, 1, nbLevels - 1 ) );

            // Convergence measures
            convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

            if ( parallel )
              convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

            multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

            int maxUsedIterations = fineModelSolid->data.rows() * fineModelSolid->data.cols();

            if ( parallel )
              maxUsedIterations += fineModelFluid->data.rows() * fineModelFluid->data.cols();

            if ( *fsiSolver == "OSM" )
              postProcessing = shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );
            else
              postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

            fineModel = shared_ptr<ImplicitMultiLevelFsiSolver> ( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

            solid.reset();
            multiLevelFluidSolver.reset();
            multiLevelSolidSolver.reset();
            convergenceMeasures.reset();
            multiLevelFsiSolver.reset();
            postProcessing.reset();

            for ( int level = 0; level < nbLevels - 1; level++ )
            {
              if ( level == 0 )
              {
                minIter = 2;
                tol = 1.0e-6;
                N = couplingGridSize / 10;
              }

              if ( level == 1 && nbLevels == 2 )
              {
                minIter = 1;
                tol = 1.0e-5;
                N = couplingGridSize;
              }

              if ( level == 1 && nbLevels == 3 )
              {
                minIter = 2;
                tol = 1.0e-6;
                N = 4 * couplingGridSize / 10;
              }

              if ( level == 2 )
              {
                minIter = 1;
                tol = 1.0e-5;
                N = couplingGridSize;
              }

              fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
              solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );
              multiLevelFluidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( fluid, fineModelFluid, 0, level ) );
              multiLevelSolidSolver = shared_ptr<MultiLevelSolver> ( new MultiLevelSolver( solid, fineModelSolid, 1, level ) );

              // Convergence measures
              convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

              convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
              convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

              if ( parallel )
                convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

              multiLevelFsiSolver = shared_ptr<MultiLevelFsiSolver> ( new MultiLevelFsiSolver( multiLevelFluidSolver, multiLevelSolidSolver, convergenceMeasures, parallel, extrapolation ) );

              int maxUsedIterations = solid->data.rows() * solid->data.cols();

              if ( parallel )
                maxUsedIterations += fluid->data.rows() * fluid->data.cols();

              if ( *fsiSolver == "OSM" )
                postProcessing = shared_ptr<PostProcessing> ( new BroydenPostProcessing( multiLevelFsiSolver, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );
              else
                postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( multiLevelFsiSolver, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

              shared_ptr<ImplicitMultiLevelFsiSolver> implicitMultiLevelFsiSolver( new ImplicitMultiLevelFsiSolver( multiLevelFsiSolver, postProcessing ) );

              models->push_back( implicitMultiLevelFsiSolver );

              solid.reset();
              multiLevelFluidSolver.reset();
              multiLevelSolidSolver.reset();
              convergenceMeasures.reset();
              multiLevelFsiSolver.reset();
              postProcessing.reset();
              implicitMultiLevelFsiSolver.reset();
            }

            models->push_back( fineModel );

            assert( static_cast<int>( models->size() ) == nbLevels );

            if ( *fsiSolver == "MM" || *fsiSolver == "OSM" || *fsiSolver == "ASM" || *fsiSolver == "ASM-ILS" )
            {
              for ( int level = 0; level < nbLevels - 1; level++ )
              {
                shared_ptr<ImplicitMultiLevelFsiSolver> fineModel;
                shared_ptr<ImplicitMultiLevelFsiSolver> coarseModel;

                fineModel = models->at( level + 1 );

                if ( level == 0 )
                  coarseModel = models->at( level );

                if ( level > 0 )
                  coarseModel = solvers->at( level - 1 );

                shared_ptr<SpaceMapping> spaceMapping;

                if ( *fsiSolver == "MM" )
                  spaceMapping = shared_ptr<SpaceMapping> ( new ManifoldMapping( fineModel, coarseModel, maxIter, singularityLimit, nbReuse, reuseInformationStartingFromTimeIndex ) );

                if ( *fsiSolver == "OSM" )
                  spaceMapping = shared_ptr<SpaceMapping> ( new OutputSpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

                if ( *fsiSolver == "ASM" )
                  spaceMapping = shared_ptr<SpaceMapping> ( new AggressiveSpaceMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex ) );

                if ( *fsiSolver == "ASM-ILS" )
                  spaceMapping = shared_ptr<SpaceMapping> ( new ASMILS( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit ) );

                shared_ptr<SpaceMappingSolver > spaceMappingSolver( new SpaceMappingSolver( fineModel, coarseModel, spaceMapping ) );

                solvers->push_back( spaceMappingSolver );

                spaceMapping.reset();
                spaceMappingSolver.reset();
              }

              assert( static_cast<int>( solvers->size() ) == nbLevels - 1 );
            }
          }

          if ( *fsiSolver == "IQN-ILS" || *fsiSolver == "IQN-Broyden" )
          {
            assert( !fsi );

            tol = 1.0e-5;
            N = couplingGridSize;
            minIter = 1;

            fluid = shared_ptr<TubeFlowFluidSolver> ( new TubeFlowFluidSolver( a0, u0, p0, dt, cmk, N, L, T, rho ) );
            solid = shared_ptr<TubeFlowSolidSolver> ( new TubeFlowSolidSolver( a0, cmk, p0, rho, L, N ) );

            shared_ptr<MultiLevelSolver> fluidSolver( new MultiLevelSolver( fluid, fluid, 0, 0 ) );
            shared_ptr<MultiLevelSolver> solidSolver( new MultiLevelSolver( solid, fluid, 1, 0 ) );

            // Convergence measures
            convergenceMeasures = shared_ptr<list<shared_ptr<ConvergenceMeasure> > >( new list<shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new MinIterationConvergenceMeasure( 0, minIter ) ) );
            convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 0, tol ) ) );

            if ( parallel )
              convergenceMeasures->push_back( shared_ptr<ConvergenceMeasure> ( new RelativeConvergenceMeasure( 1, tol ) ) );

            int maxUsedIterations = solid->data.rows() * solid->data.cols();

            if ( parallel )
              maxUsedIterations += fluid->data.rows() * fluid->data.cols();

            fsi = shared_ptr<MultiLevelFsiSolver>( new MultiLevelFsiSolver( fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation ) );

            if ( *fsiSolver == "IQN-ILS" )
              postProcessing = shared_ptr<PostProcessing> ( new IQNILSPostProcessing( fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling ) );

            if ( *fsiSolver == "IQN-Broyden" )
              postProcessing = shared_ptr<PostProcessing> ( new BroydenPostProcessing( fsi, initialRelaxation, maxIter, maxUsedIterations, nbReuse, reuseInformationStartingFromTimeIndex ) );

            solver = shared_ptr<Solver>( new ImplicitMultiLevelFsiSolver( fsi, postProcessing ) );
          }

          if ( *fsiSolver == "MM" || *fsiSolver == "OSM" || *fsiSolver == "ASM" || *fsiSolver == "ASM-ILS" )
            solver = shared_ptr<Solver> ( new MultiLevelSpaceMappingSolver( solvers, models ) );

          if ( *fsiSolver == "ML-IQN-ILS" )
            solver = shared_ptr<Solver> ( new MLIQNILSSolver( models ) );

          try
          {
            solver->run();
          }
          catch ( ... )
          {
            continue;
          }
          ofstream logFile( label + ".log" );
          assert( logFile.is_open() );

          logFile << "label = " << label << endl;
          logFile << "nbLevels = " << nbLevels << endl;
          logFile << "solver = " << *fsiSolver << endl;
          logFile << "nbReuse = " << nbReuse << endl;
          logFile << "parallel = " << parallel << endl;
          logFile << "elapsed time = " << solver->timeElapsed() << endl;
          logFile << "nbTimeSteps = " << fluid->timeIndex << endl;

          if ( *fsiSolver == "IQN-ILS" || *fsiSolver == "IQN-Broyden" )
            logFile << "nbIter = " << fsi->nbIter << endl;

          if ( *fsiSolver == "MM" || *fsiSolver == "OSM" || *fsiSolver == "ML-IQN-ILS" || *fsiSolver == "ASM" || *fsiSolver == "ASM-ILS" )
          {
            int level = 0;

            for ( std::deque<shared_ptr<ImplicitMultiLevelFsiSolver> >::iterator it = models->begin(); it != models->end(); ++it )
            {
              shared_ptr<ImplicitMultiLevelFsiSolver> model = *it;
              logFile << "level " << level << " nbIter = " << model->fsi->nbIter << endl;

              level++;
            }
          }

          logFile.close();
        }
      }
    }
  }

  return 0;
}
