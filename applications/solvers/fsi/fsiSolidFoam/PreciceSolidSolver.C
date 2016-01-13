
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PreciceSolidSolver.H"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixRowMajor;

PreciceSolidSolver::PreciceSolidSolver( shared_ptr<foamSolidSolver> solver )
    :
    solver( solver ),
    precice( shared_ptr<precice::SolverInterface> ( new precice::SolverInterface( "Structure_Solver", Pstream::myProcNo(), Pstream::nProcs() ) ) ),
    idsReadPositions(),
    idsWritePositions()
{
    assert( solver );

    init();
}

PreciceSolidSolver::~PreciceSolidSolver()
{
    precice->finalize();
}

void PreciceSolidSolver::init()
{
    string filename = static_cast<std::string>( solver->args->rootPath() ) + "/" + static_cast<std::string>( solver->args->globalCaseName() ) + "/constant/preCICE.xml";
    precice->configure( filename );

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    assert( precice->getDimensions() == solver->mesh.nGeometricD() );

    setReadPositions();
    setWritePositions();

    precice->initialize();
}

void PreciceSolidSolver::readData( matrix & data )
{
    // Read displacements from preCICE

    int meshId = precice->getMeshID( "Structure_CellCenters" );
    int dataId = precice->getDataID( "Stresses", meshId );

    matrixRowMajor dataRowMajor( idsReadPositions.rows(), precice->getDimensions() );
    dataRowMajor.setZero();

    if ( data.rows() > 0 )
        precice->readBlockVectorData( dataId, idsReadPositions.rows(), idsReadPositions.data(), dataRowMajor.data() );

    data = dataRowMajor.cast<scalar>();
}

void PreciceSolidSolver::run()
{
    matrix outputOld( solver->getInterfaceSizeLocal(), precice->getDimensions() );
    matrix input( solver->getInterfaceSizeLocal(), precice->getDimensions() );
    matrix output;
    outputOld.setZero();

    while ( solver->isRunning() )
    {
        solver->initTimeStep();

        int iter = 0;

        while ( precice->isCouplingOngoing() )
        {
            Info << endl << "Time = " << solver->runTime->timeName() << ", iteration = " << iter + 1 << endl;

            readData( input );

            if ( precice->isActionRequired( precice::constants::actionReadIterationCheckpoint() ) )
                precice->fulfilledAction( precice::constants::actionReadIterationCheckpoint() );

            solver->setTractionLocal( input );
            solver->solve();
            solver->getDisplacementLocal( output );

            writeData( output - outputOld );

            if ( precice->isActionRequired( precice::constants::actionWriteIterationCheckpoint() ) )
                precice->fulfilledAction( precice::constants::actionWriteIterationCheckpoint() );

            precice->advance( solver->runTime->deltaT().value() );

            iter++;

            if ( precice->isTimestepComplete() )
                break;
        }

        solver->finalizeTimeStep();

        outputOld = output;
    }
}

void PreciceSolidSolver::setReadPositions()
{
    // Initialize matrices
    matrix readPositionsColumnMajor;
    matrixRowMajor readPositions;

    // Retrieve positions from solver solver
    solver->getReadPositionsLocal( readPositionsColumnMajor );

    // Store the positions in row-major for preCICE. Eigen uses column major by default.
    readPositions = readPositionsColumnMajor.cast<double>();

    assert( readPositions.cols() == precice->getDimensions() );

    // Get the mesh id
    int meshId = precice->getMeshID( "Structure_CellCenters" );

    // Send the positions to preCICE

    idsReadPositions.resize( readPositions.rows() );
    idsReadPositions.setZero();

    if ( readPositions.rows() > 0 )
        precice->setMeshVertices( meshId, readPositions.rows(), readPositions.data(), idsReadPositions.data() );
}

void PreciceSolidSolver::setWritePositions()
{
    // Initialize matrices
    matrix writePositionsColumnMajor;
    matrixRowMajor writePositions;

    // Retrieve positions from solver solver
    solver->getWritePositionsLocal( writePositionsColumnMajor );

    // Store the positions in row-major for preCICE. Eigen uses column major by default.
    writePositions = writePositionsColumnMajor.cast<double>();

    assert( writePositions.cols() == precice->getDimensions() );

    // Get the mesh id
    int meshId = precice->getMeshID( "Structure_Nodes" );

    // Send the positions to preCICE

    idsWritePositions.resize( writePositions.rows() );
    idsWritePositions.setZero();

    if ( writePositions.rows() > 0 )
        precice->setMeshVertices( meshId, writePositions.rows(), writePositions.data(), idsWritePositions.data() );
}

void PreciceSolidSolver::writeData( const matrix & data )
{
    // Send forces to preCICE
    matrixRowMajor dataRowMajor = data.cast<double>();

    int meshId = precice->getMeshID( "Structure_Nodes" );
    int dataId = precice->getDataID( "Displacements", meshId );

    assert( data.cols() == precice->getDimensions() );

    if ( dataRowMajor.rows() > 0 )
        precice->writeBlockVectorData( dataId, dataRowMajor.rows(), idsWritePositions.data(), dataRowMajor.data() );
}
