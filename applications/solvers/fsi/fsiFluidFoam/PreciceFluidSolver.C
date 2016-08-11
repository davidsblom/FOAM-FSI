
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
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    int meshId = precice->getMeshID( "Fluid_Acoustics" );
    int dataIdPressure = precice->getDataID( "Acoustics_Pressure", meshId );
    int dataIdDensity = precice->getDataID( "Acoustics_Density", meshId );
    int dataIdVelocityX = precice->getDataID( "Acoustics_Velocity_X", meshId );
    int dataIdVelocityY = precice->getDataID( "Acoustics_Velocity_Y", meshId );

    int dataIdVelocityZ = -1;

    if ( precice->hasData( "Acoustics_Velocity_Z", meshId ) )
        dataIdVelocityZ = precice->getDataID( "Acoustics_Velocity_Z", meshId );

    matrix velocityData( idsWritePositionsAcoustics.size(), precice->getDimensions() );

    if ( velocityData.rows() > 0 && precice->hasData( "Acoustics_Velocity_X", meshId ) )
    {
        matrixRowMajor dataVelocityX( velocityData.rows(), 1 ), dataVelocityY( velocityData.rows(), 1 );
        precice->readBlockScalarData( dataIdVelocityX, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataVelocityX.data() );
        precice->readBlockScalarData( dataIdVelocityY, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataVelocityY.data() );
        velocityData.col( 0 ) = dataVelocityX.cast<scalar>();
        velocityData.col( 1 ) = dataVelocityY.cast<scalar>();

        if ( precice->hasData( "Acoustics_Velocity_Z", meshId ) )
        {
            matrixRowMajor dataVelocityZ( velocityData.rows(), 1 );
            precice->readBlockScalarData( dataIdVelocityZ, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), dataVelocityZ.data() );
            velocityData.col( 2 ) = dataVelocityZ.cast<scalar>();
        }
    }

    matrixRowMajor pressureDataRowMajor( idsWritePositionsAcoustics.size(), 1 );

    if ( pressureDataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Pressure", meshId ) )
        precice->readBlockScalarData( dataIdPressure, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), pressureDataRowMajor.data() );

    matrixRowMajor densityDataRowMajor( idsWritePositionsAcoustics.size(), 1 );

    if ( densityDataRowMajor.rows() > 0 && precice->hasData( "Acoustics_Density", meshId ) )
        precice->readBlockScalarData( dataIdDensity, idsWritePositionsAcoustics.rows(), idsWritePositionsAcoustics.data(), densityDataRowMajor.data() );

    const matrix densityData = densityDataRowMajor.cast<scalar>();
    const matrix pressureData = pressureDataRowMajor.cast<scalar>();

    // Perfect gas law: p = rho * R * T, T = p / ( R * rho )
    const double R = 296.0; // N m / kmol K
    assert( pressureData.rows() == densityData.rows() );
    matrix temperatureData = pressureData.array() / (R * densityData.array() + SMALL);

    solver->setAcousticsData( pressureData, temperatureData, velocityData );
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
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
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
    Info << "Fluid-Acoustics interface: " << sum( interfaceSize ) << " points" << endl;

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
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    const int meshId = precice->getMeshID( "Fluid_Acoustics" );

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

void PreciceFluidSolver::writeAcousticsDensityGradient()
{
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    const int meshId = precice->getMeshID( "Fluid_Acoustics" );

    if ( !precice->hasData( "Acoustics_Density_Gradient", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Density_Gradient", meshId );

    matrix data;
    solver->getAcousticsDensityGradientLocal( data );

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsPressure()
{
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    const int meshId = precice->getMeshID( "Fluid_Acoustics" );

    if ( !precice->hasData( "Acoustics_Density", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Density", meshId );

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
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    const int meshId = precice->getMeshID( "Fluid_Acoustics" );

    if ( !precice->hasData( "Acoustics_Pressure_Gradient", meshId ) )
        return;

    const int dataId = precice->getDataID( "Acoustics_Pressure_Gradient", meshId );

    matrix data;
    solver->getAcousticsPressureGradientLocal( data );

    if ( data.rows() > 0 )
    {
        matrixRowMajor dataRowMajor = data.cast<double>();
        precice->writeBlockScalarData( dataId, data.rows(), idsWritePositionsAcoustics.data(), dataRowMajor.data() );
    }
}

void PreciceFluidSolver::writeAcousticsVelocity()
{
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    const int meshId = precice->getMeshID( "Fluid_Acoustics" );

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
    if ( !precice->hasMesh( "Fluid_Acoustics" ) )
        return;

    const int meshId = precice->getMeshID( "Fluid_Acoustics" );

    matrix data;
    solver->getAcousticsVelocityGradientLocal( data );

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

void PreciceFluidSolver::writeDataAcoustics()
{
    // writeAcousticsPressure();
    // writeAcousticsDensity();
    // writeAcousticsVelocity();
    writeAcousticsDensityGradient();
    writeAcousticsPressureGradient();
    writeAcousticsVelocityGradient();
}
