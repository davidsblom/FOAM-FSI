
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
    idsWritePositions(),
    totalRunTime( 0 ),
    totalNbIterations( 0 )
{
    assert( solver );

    init();
}

PreciceFluidSolver::~PreciceFluidSolver()
{
    assert( !precice->isCouplingOngoing() );
    precice->finalize();
}

void PreciceFluidSolver::init()
{
    string filename = static_cast<std::string>( solver->args->rootPath() ) + "/" + static_cast<std::string>( solver->args->globalCaseName() ) + "/constant/preCICE.xml";
    precice->configure( filename );

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    assert( precice->getDimensions() == solver->mesh.nGeometricD() );

    setReadPositions();
    setWritePositions();
    setWritePositionsAcoustics();

    precice->initialize();

    writeDataAcoustics();

    precice->fulfilledAction( precice::constants::actionWriteInitialData() );

    precice->initializeData();
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
            std::clock_t t = std::clock();

            Info << endl << "Time = " << solver->runTime->timeName() << ", iteration = " << iter + 1 << endl;

            matrix input( solver->getInterfaceSizeLocal(), precice->getDimensions() ), output;

            input.setZero();

            readData( input );

            if ( precice->isActionRequired( precice::constants::actionReadIterationCheckpoint() ) )
                precice->fulfilledAction( precice::constants::actionReadIterationCheckpoint() );

            solver->setDisplacementLocal( input );
            solver->moveMesh();
            solver->solve();
            solver->getTractionLocal( output );

            writeData( output );
            writeDataAcoustics();

            if ( precice->isActionRequired( precice::constants::actionWriteIterationCheckpoint() ) )
                precice->fulfilledAction( precice::constants::actionWriteIterationCheckpoint() );

            t = std::clock() - t;
            double runTime = static_cast<float>(t) / CLOCKS_PER_SEC;
            totalRunTime += runTime;
            totalNbIterations++;
            Info << "runtime = " << runTime << " s" << endl;
            Info << "average runtime = " << totalRunTime / totalNbIterations << " s" << endl;
            Info << "total runtime = " << totalRunTime << " s" << endl;

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

    assert( readPositionsColumnMajor.cols() == precice->getDimensions() );

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

void PreciceFluidSolver::setWritePositionsAcoustics()
{
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    // Initialize matrices
    matrix writePositionsColumnMajor;
    matrixRowMajor writePositions;

    // Retrieve positions from fluid solver
    solver->getWritePositionsLocalAcoustics( writePositionsColumnMajor );

    assert( writePositionsColumnMajor.cols() == precice->getDimensions() );

    // Store the positions in row-major order for preCICE. Eigen uses column major by default
    writePositions = writePositionsColumnMajor;

    // Get the mesh id
    int meshId = precice->getMeshID( "Fluid_Acoustics" );

    // Send the write positions to preCICE
    idsWritePositionsAcoustics.resize( writePositions.rows() );
    idsWritePositionsAcoustics.setZero();

    if ( writePositions.rows() > 0 )
        precice->setMeshVertices( meshId, writePositions.rows(), writePositions.data(), idsWritePositionsAcoustics.data() );
}

void PreciceFluidSolver::setWritePositions()
{
    // Initialize matrices
    matrix writePositionsColumnMajor;
    matrixRowMajor writePositions;

    // Retrieve positions from solver solver
    solver->getWritePositionsLocal( writePositionsColumnMajor );

    assert( writePositionsColumnMajor.cols() == precice->getDimensions() );

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

    assert( data.cols() == precice->getDimensions() );

    if ( dataRowMajor.rows() > 0 )
        precice->writeBlockVectorData( dataId, dataRowMajor.rows(), idsWritePositions.data(), dataRowMajor.data() );
}

void PreciceFluidSolver::writeDataAcoustics()
{
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    // Initialize variables
    int meshId, dataIdDensity, dataIdPressure, dataIdVelocityX, dataIdVelocityY, dataIdVelocityZ;
    matrix dataDensity, dataVelocity, dataPressure;

    // Send data to preCICE

    meshId = precice->getMeshID( "Fluid_Acoustics" );
    dataIdDensity = precice->getDataID( "Acoustics_Density", meshId );
    dataIdVelocityX = precice->getDataID( "Acoustics_Velocity_X", meshId );
    dataIdVelocityY = precice->getDataID( "Acoustics_Velocity_Y", meshId );
    dataIdVelocityZ = -1;
    dataIdPressure = -1;

    if ( precice->hasData( "Acoustics_Velocity_Z", meshId ) )
        dataIdVelocityZ = precice->getDataID( "Acoustics_Velocity_Z", meshId );

    if ( precice->hasData( "Acoustics_Pressure", meshId ) )
    {
        dataIdPressure = precice->getDataID( "Acoustics_Pressure", meshId );
        solver->getAcousticsPressureLocal( dataPressure );

        assert( dataPressure.rows() == idsWritePositionsAcoustics.rows() );
    }

    solver->getAcousticsDensityLocal( dataDensity );
    solver->getAcousticsVelocityLocal( dataVelocity );

    assert( dataDensity.rows() == idsWritePositionsAcoustics.rows() );
    assert( dataVelocity.rows() == idsWritePositionsAcoustics.rows() );
    assert( dataVelocity.cols() == precice->getDimensions() );

    if ( dataDensity.rows() > 0 )
        precice->writeBlockScalarData( dataIdDensity, dataDensity.rows(), idsWritePositionsAcoustics.data(), dataDensity.data() );

    if ( dataPressure.rows() > 0 )
        precice->writeBlockScalarData( dataIdPressure, dataPressure.rows(), idsWritePositionsAcoustics.data(), dataPressure.data() );

    if ( dataVelocity.rows() > 0 )
    {
        matrix dataVelocityX, dataVelocityY, dataVelocityZ;
        dataVelocityX = dataVelocity.col( 0 );
        dataVelocityY = dataVelocity.col( 1 );
        dataVelocityZ = dataVelocity.col( 2 );

        precice->writeBlockScalarData( dataIdVelocityX, dataVelocityX.rows(), idsWritePositionsAcoustics.data(), dataVelocityX.data() );
        precice->writeBlockScalarData( dataIdVelocityY, dataVelocityY.rows(), idsWritePositionsAcoustics.data(), dataVelocityY.data() );

        if ( precice->hasData( "Acoustics_Velocity_Z", meshId ) )
        {
            assert( dataVelocity.cols() == 3 );
            assert( precice->getDimensions() == 3 );
            precice->writeBlockScalarData( dataIdVelocityZ, dataVelocityZ.rows(), idsWritePositionsAcoustics.data(), dataVelocityZ.data() );
        }
    }
}
