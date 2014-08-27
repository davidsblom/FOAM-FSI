
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PreciceFluidSolver.H"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixRowMajor;

PreciceFluidSolver::PreciceFluidSolver( shared_ptr<foamFluidSolver> solver )
  :
  solver( solver ),
  precice( shared_ptr<precice::SolverInterface> ( new precice::SolverInterface( "Fluid_Solver", Pstream::myProcNo(), Pstream::nProcs() ) ) ),
  idsReadPositions(),
  idsWritePositions()
{
  assert( solver );

  init();
}

PreciceFluidSolver::~PreciceFluidSolver()
{}

void PreciceFluidSolver::init()
{
  string filename = static_cast<std::string>( solver->args->rootPath() ) + "/" + static_cast<std::string>( solver->args->globalCaseName() ) + "/constant/preCICE.xml";
  precice->configure( filename );

  assert( precice->getDimensions() == 3 );

  setReadPositions();

  setWritePositions();

  precice->initialize();
}

void PreciceFluidSolver::readData( matrix & data )
{
  // Read displacements from preCICE

  int meshId = precice->getMeshID( "Fluid_Nodes" );
  int dataId = precice->getDataID( "Displacements", meshId );

  matrixRowMajor dataRowMajor( idsReadPositions.rows(), precice->getDimensions() );
  dataRowMajor.setZero();

  if ( data.rows() > 0 )
    precice->readBlockVectorData( dataId, idsReadPositions.rows(), idsReadPositions.data(), dataRowMajor.data() );

  data = dataRowMajor;
}

void PreciceFluidSolver::run()
{
  while ( solver->isRunning() )
  {
    solver->initTimeStep();

    int iter = 0;

    while ( precice->isCouplingOngoing() )
    {
      Info << endl << "Time = " << solver->runTime->timeName() << ", iteration = " << iter + 1 << endl;

      matrix input( solver->getInterfacePointsSizeLocal(), precice->getDimensions() ), output;

      input.setZero();

      if ( precice->isReadDataAvailable() )
        readData( input );

      if ( precice->isActionRequired( precice::constants::actionReadIterationCheckpoint() ) )
        precice->fulfilledAction( precice::constants::actionReadIterationCheckpoint() );

      solver->setDisplacementLocal( input );
      solver->solve();
      solver->getTractionLocal( output );

      if ( precice->isWriteDataRequired( solver->runTime->deltaT().value() ) )
        writeData( output );

      if ( precice->isActionRequired( precice::constants::actionWriteIterationCheckpoint() ) )
        precice->fulfilledAction( precice::constants::actionWriteIterationCheckpoint() );

      precice->advance( solver->runTime->deltaT().value() );

      iter++;

      if ( precice->isTimestepComplete() )
        break;
    }

    solver->finalizeTimeStep();
  }
}

void PreciceFluidSolver::setReadPositions()
{
  // Initialize matrices
  matrix readPositionsColumnMajor;
  matrixRowMajor readPositions;

  // Retrieve positions from solver solver
  solver->getReadPositionsLocal( readPositionsColumnMajor );

  // Store the positions in row-major for preCICE. Eigen uses column major by default.
  readPositions = readPositionsColumnMajor;

  // Get the mesh id
  int meshId = precice->getMeshID( "Fluid_Nodes" );

  // Send the positions to preCICE

  idsReadPositions.resize( readPositions.rows() );
  idsReadPositions.setZero();

  if ( readPositions.rows() > 0 )
    precice->setMeshVertices( meshId, readPositions.rows(), readPositions.data(), idsReadPositions.data() );
}

void PreciceFluidSolver::setWritePositions()
{
  // Initialize matrices
  matrix writePositionsColumnMajor;
  matrixRowMajor writePositions;

  // Retrieve positions from solver solver
  solver->getWritePositionsLocal( writePositionsColumnMajor );

  // Store the positions in row-major for preCICE. Eigen uses column major by default
  writePositions = writePositionsColumnMajor;

  // Get the mesh id
  int meshId = precice->getMeshID( "Fluid_CellCenters" );

  // Send the write positions to preCICE

  idsWritePositions.resize( writePositions.rows() );
  idsWritePositions.setZero();

  if ( writePositions.rows() > 0 )
    precice->setMeshVertices( meshId, writePositions.rows(), writePositions.data(), idsWritePositions.data() );
}

void PreciceFluidSolver::writeData( const matrix & data )
{
  // Send forces to preCICE
  matrixRowMajor dataRowMajor = data;

  int meshId = precice->getMeshID( "Fluid_CellCenters" );
  int dataId = precice->getDataID( "Stresses", meshId );

  if ( dataRowMajor.rows() > 0 )
    precice->writeBlockVectorData( dataId, dataRowMajor.rows(), idsWritePositions.data(), dataRowMajor.data() );
}
