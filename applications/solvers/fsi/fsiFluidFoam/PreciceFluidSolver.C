
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PreciceFluidSolver.H"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixRowMajor;

PreciceFluidSolver::PreciceFluidSolver(
    shared_ptr<foamFluidSolver> solver,
    const std::string & participant
    )
    :
    solver( solver ),
    precice( shared_ptr<precice::SolverInterface> ( new precice::SolverInterface( participant, Pstream::myProcNo(), Pstream::nProcs() ) ) ),
    idsReadPositions(),
    idsWritePositions(),
    FLUID_ACOUSTICS_WRITE( "Fluid_Acoustics" ),
    FLUID_ACOUSTICS_READ( "Fluid_Acoustics" ),
    participant( participant )
{
    assert( solver );
    assert( participant == "Fluid_Solver" || participant == "Fluid_Solver_A" || participant == "Fluid_Solver_B" );

    if ( participant == "Fluid_Solver_A" )
        FLUID_ACOUSTICS_WRITE = "Fluid_Acoustics_A";

    if ( participant == "Fluid_Solver_B" )
        FLUID_ACOUSTICS_WRITE = "Fluid_Acoustics_B";

    FLUID_ACOUSTICS_READ = FLUID_ACOUSTICS_WRITE;

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

    try
    {
        precice->configure( filename );
    }
    catch ( std::string e )
    {
        Info << "Error: initialization of precice failed" << nl << e << endl;
        abort( FatalError );
    }
    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    assert( precice->getDimensions() == solver->mesh.nGeometricD() );

    setReadPositions();
    setWritePositions();
    setWritePositionsAcoustics();

    writeDataAcoustics();

    precice->initialize();

    precice->fulfilledAction( precice::constants::actionWriteInitialData() );

    precice->initializeData();
}

void PreciceFluidSolver::readData( matrix & data )
{
    if ( !precice->hasMesh( "Fluid_Nodes" ) )
        return;

    // Read displacements from preCICE

    int meshId = precice->getMeshID( "Fluid_Nodes" );
    int dataId = precice->getDataID( "Displacements", meshId );

    matrixRowMajor dataRowMajor( idsReadPositions.rows(), precice->getDimensions() );
    dataRowMajor.setZero();

    if ( dataRowMajor.rows() > 0 )
        precice->readBlockVectorData( dataId, idsReadPositions.rows(), idsReadPositions.data(), dataRowMajor.data() );

    data = dataRowMajor.cast<scalar>();
}

void PreciceFluidSolver::readDataAcoustics()
{
    if ( participant == "Fluid_Solver_A" )
    {
        readDataAcousticsTemperature();
        readDataAcousticsPressure();
        readDataAcousticsVelocity();
    }

    if ( participant == "Fluid_Solver_B" )
    {
        readDataAcousticsTemperatureGradient();
        readDataAcousticsPressureGradient();
        readDataAcousticsVelocityGradient();
    }
}

void PreciceFluidSolver::writeDataAcoustics()
{
    if ( participant == "Fluid_Solver_B" )
    {
        writeAcousticsDensity();
        writeAcousticsPressure();
        writeAcousticsVelocity();
        writeAcousticsTemperature();
    }

    if ( participant == "Fluid_Solver_A" )
    {
        writeAcousticsPressureGradient();
        writeAcousticsVelocityGradient();
        writeAcousticsTemperatureGradient();
    }
}

void PreciceFluidSolver::readDataAcousticsTemperature()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_READ ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_READ );

    if ( !precice->hasData( "Acoustics_Temperature", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Temperature", meshId );

    matrixRowMajor dataRowMajor( idsWritePositionsAcoustics.size(), 1 );

    if ( dataRowMajor.rows() > 0 )
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );

    matrix data = dataRowMajor.cast<scalar>();
    solver->setAcousticsTemperature( data );
}

void PreciceFluidSolver::readDataAcousticsTemperatureGradient()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_READ ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_READ );

    if ( !precice->hasData( "Acoustics_Temperature_Gradient", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Temperature_Gradient", meshId );

    matrixRowMajor dataRowMajor( idsWritePositionsAcoustics.size(), 1 );

    if ( dataRowMajor.rows() > 0 )
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );

    matrix data = dataRowMajor.cast<scalar>();
    solver->setAcousticsTemperatureGradient( data );
}

void PreciceFluidSolver::readDataAcousticsPressure()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_READ ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_READ );

    if ( !precice->hasData( "Acoustics_Pressure", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Pressure", meshId );

    matrixRowMajor dataRowMajor( idsWritePositionsAcoustics.size(), 1 );

    if ( dataRowMajor.rows() > 0 )
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );

    matrix data = dataRowMajor.cast<scalar>();
    solver->setAcousticsPressure( data );
}

void PreciceFluidSolver::readDataAcousticsPressureGradient()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_READ ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_READ );

    if ( !precice->hasData( "Acoustics_Pressure_Gradient", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Pressure_Gradient", meshId );

    matrixRowMajor dataRowMajor( idsWritePositionsAcoustics.size(), 1 );

    if ( dataRowMajor.rows() > 0 )
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );

    matrix data = dataRowMajor.cast<scalar>();
    solver->setAcousticsPressureGradient( data );
}

void PreciceFluidSolver::readDataAcousticsVelocity()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_READ ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_READ );

    matrixRowMajor dataRowMajor( idsWritePositionsAcoustics.size(), precice->getDimensions() );

    if ( dataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Velocity_X", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_X", meshId );
        matrixRowMajor dataX( dataRowMajor.rows(), 1 );
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataX.data() );
        dataRowMajor.col( 0 ) = dataX;
    }

    if ( dataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Velocity_Y", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Y", meshId );
        matrixRowMajor dataY( dataRowMajor.rows(), 1 );
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataY.data() );
        dataRowMajor.col( 1 ) = dataY;
    }

    if ( dataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Velocity_Z", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Z", meshId );
        matrixRowMajor dataZ( dataRowMajor.rows(), 1 );
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataZ.data() );
        dataRowMajor.col( 2 ) = dataZ;
    }

    if ( precice->hasData( "Acoustics_Velocity_X", meshId ) && precice->hasData( "Acoustics_Velocity_Y", meshId ) )
    {
        matrix data = dataRowMajor.cast<scalar>();
        solver->setAcousticsVelocity( data );
    }
}

void PreciceFluidSolver::readDataAcousticsVelocityGradient()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_READ ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_READ );

    matrixRowMajor dataRowMajor( idsWritePositionsAcoustics.size(), precice->getDimensions() );

    if ( dataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Velocity_X_Gradient", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_X_Gradient", meshId );
        matrixRowMajor dataX( dataRowMajor.rows(), 1 );
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataX.data() );
        dataRowMajor.col( 0 ) = dataX;
    }

    if ( dataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Velocity_Y_Gradient", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Y_Gradient", meshId );
        matrixRowMajor dataY( dataRowMajor.rows(), 1 );
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataY.data() );
        dataRowMajor.col( 1 ) = dataY;
    }

    if ( dataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Velocity_Z_Gradient", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Z_Gradient", meshId );
        matrixRowMajor dataZ( dataRowMajor.rows(), 1 );
        precice->readBlockScalarData( dataId, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataZ.data() );
        dataRowMajor.col( 2 ) = dataZ;
    }

    if ( precice->hasData( "Acoustics_Velocity_X_Gradient", meshId ) && precice->hasData( "Acoustics_Velocity_Y_Gradient", meshId ) )
    {
        matrix data = dataRowMajor.cast<scalar>();
        solver->setAcousticsVelocityGradient( data );
    }
}

void PreciceFluidSolver::run()
{
    matrix input, inputOld, output;

    while ( solver->isRunning() )
    {
        solver->initTimeStep();

        int iter = 0;

        while ( precice->isCouplingOngoing() )
        {
            Info << endl << "Time = " << solver->runTime->timeName() << ", iteration = " << iter + 1 << endl;

            readData( input );
            readDataAcoustics();

            if ( precice->isActionRequired( precice::constants::actionReadIterationCheckpoint() ) )
                precice->fulfilledAction( precice::constants::actionReadIterationCheckpoint() );

            if ( precice->hasMesh( "Fluid_Nodes" ) )
            {
                if ( input.cols() == inputOld.cols() )
                    solver->setDisplacementLocal( input + inputOld );
                else
                    solver->setDisplacementLocal( input );

                solver->moveMesh();
            }

            solver->solve();

            if ( precice->hasMesh( "Fluid_CellCenters" ) )
                solver->getTractionLocal( output );

            writeData( output );
            writeDataAcoustics();

            if ( precice->isActionRequired( precice::constants::actionWriteIterationCheckpoint() ) )
                precice->fulfilledAction( precice::constants::actionWriteIterationCheckpoint() );

            precice->advance( solver->runTime->deltaT().value() );

            iter++;

            if ( precice->isTimestepComplete() )
                break;
        }

        solver->finalizeTimeStep();

        if ( input.cols() == inputOld.cols() )
            inputOld += input;
        else
            inputOld = input;

        if ( not precice->isCouplingOngoing() )
            break;
    }
}

void PreciceFluidSolver::setReadPositions()
{
    if ( !precice->hasMesh( "Fluid_Nodes" ) )
        return;

    // Initialize matrices
    matrix readPositionsColumnMajor;
    matrixRowMajor readPositions;

    // Retrieve positions from solver solver
    solver->getReadPositionsLocal( readPositionsColumnMajor );

    assert( readPositionsColumnMajor.cols() == precice->getDimensions() );

    // Store the positions in row-major for preCICE. Eigen uses column major by default.
    readPositions = readPositionsColumnMajor.cast<double>();

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
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    // Initialize matrices
    matrix writePositionsColumnMajor;
    matrixRowMajor writePositions;

    // Retrieve positions from fluid solver
    solver->getWritePositionsLocalAcoustics( writePositionsColumnMajor );

    assert( writePositionsColumnMajor.cols() == precice->getDimensions() );

    // Store the positions in row-major order for preCICE. Eigen uses column major by default
    writePositions = writePositionsColumnMajor.cast<double>();

    labelList interfaceSize( Pstream::nProcs(), 0 );
    interfaceSize[Pstream::myProcNo()] = writePositions.rows();
    reduce( interfaceSize, sumOp<labelList>() );
    Info << "Fluid-Acoustics-Write mesh size: " << sum( interfaceSize ) << " points" << endl;

    // Get the mesh id
    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    // Send the write positions to preCICE
    idsWritePositionsAcoustics.resize( writePositions.rows() );
    idsWritePositionsAcoustics.setZero();

    if ( writePositions.rows() > 0 )
        precice->setMeshVertices( meshId, writePositions.rows(), writePositions.data(), idsWritePositionsAcoustics.data() );
}

void PreciceFluidSolver::setWritePositions()
{
    if ( !precice->hasMesh( "Fluid_CellCenters" ) )
        return;

    // Initialize matrices
    matrix writePositionsColumnMajor;
    matrixRowMajor writePositions;

    // Retrieve positions from solver solver
    solver->getWritePositionsLocal( writePositionsColumnMajor );

    assert( writePositionsColumnMajor.cols() == precice->getDimensions() );

    // Store the positions in row-major for preCICE. Eigen uses column major by default
    writePositions = writePositionsColumnMajor.cast<double>();

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
    if ( !precice->hasMesh( "Fluid_CellCenters" ) )
        return;

    // Send forces to preCICE
    matrixRowMajor dataRowMajor = data.cast<double>();

    int meshId = precice->getMeshID( "Fluid_CellCenters" );
    int dataId = precice->getDataID( "Stresses", meshId );

    assert( data.cols() == precice->getDimensions() );

    if ( dataRowMajor.rows() > 0 )
        precice->writeBlockVectorData( dataId, dataRowMajor.rows(), idsWritePositions.data(), dataRowMajor.data() );
}

void PreciceFluidSolver::writeAcousticsDensity()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    if ( !precice->hasData( "Acoustics_Density", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Density", meshId );

    matrix data;
    solver->getAcousticsDensityLocal( data );

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsPressure()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    if ( !precice->hasData( "Acoustics_Pressure", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Pressure", meshId );

    matrix data;
    solver->getAcousticsPressureLocal( data );

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsPressureGradient()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    if ( !precice->hasData( "Acoustics_Pressure_Gradient", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Pressure_Gradient", meshId );

    matrix data;
    solver->getAcousticsPressureGradientLocal( data );
    data *= -1.0;

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsTemperature()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    if ( !precice->hasData( "Acoustics_Temperature", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Temperature", meshId );

    matrix data;
    solver->getAcousticsTemperatureLocal( data );

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsTemperatureGradient()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    if ( !precice->hasData( "Acoustics_Temperature_Gradient", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Temperature_Gradient", meshId );

    matrix data;
    solver->getAcousticsTemperatureGradientLocal( data );
    data *= -1.0;

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsVelocity()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    matrix data;
    solver->getAcousticsVelocityLocal( data );

    if ( data.rows() > 0 && precice->hasData( "Acoustics_Velocity_X", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_X", meshId );
        matrixRowMajor dataRowMajor = data.col( 0 ).cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }

    if ( data.rows() > 0 && precice->hasData( "Acoustics_Velocity_Y", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Y", meshId );
        matrixRowMajor dataRowMajor = data.col( 1 ).cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }

    if ( data.rows() > 0 && precice->hasData( "Acoustics_Velocity_Z", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Z", meshId );
        matrixRowMajor dataRowMajor = data.col( 2 ).cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsVelocityGradient()
{
    if ( !precice->hasMesh( FLUID_ACOUSTICS_WRITE ) )
        return;

    const int meshId = precice->getMeshID( FLUID_ACOUSTICS_WRITE );

    matrix data;
    solver->getAcousticsVelocityGradientLocal( data );
    data *= -1.0;

    if ( data.rows() > 0 && precice->hasData( "Acoustics_Velocity_X_Gradient", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_X_Gradient", meshId );
        matrixRowMajor dataRowMajor = data.col( 0 ).cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }

    if ( data.rows() > 0 && precice->hasData( "Acoustics_Velocity_Y_Gradient", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Y_Gradient", meshId );
        matrixRowMajor dataRowMajor = data.col( 1 ).cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }

    if ( data.rows() > 0 && precice->hasData( "Acoustics_Velocity_Z_Gradient", meshId ) )
    {
        const int dataId = precice->getDataID( "Acoustics_Velocity_Z_Gradient", meshId );
        matrixRowMajor dataRowMajor = data.col( 2 ).cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}
