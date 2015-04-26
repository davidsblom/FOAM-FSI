
#include <iostream>
#include <memory>
#include "fvCFD.H"
#include "SolidSolver.H"
#include "CoupledFluidSolver.H"
#include "RBFInterpolation.H"
#include "TPSFunction.H"

class MapMeshSolidSolver : public SolidSolver
{
public:

    MapMeshSolidSolver(
        string name,
        std::shared_ptr<argList> args,
        std::shared_ptr<Time> runTime
        )
        :
        SolidSolver( name, args, runTime )
    {}

    void getDisplacementLocal( matrix & displacement )
    {
        displacement.resize( getInterfaceSizeLocal(), mesh.nGeometricD() );

        int offset = 0;

        forAll( movingPatchIDs, patchI )
        {
            int size = U.boundaryField()[movingPatchIDs[patchI]].size();

            for ( int i = 0; i < size; i++ )
                for ( int j = 0; j < displacement.cols(); j++ )
                    displacement( i + offset, j ) = U.boundaryField()[movingPatchIDs[patchI]][i][j];

            offset += size;
        }
    }

    virtual void getWritePositions( matrix & writePositions )
    {
        vectorField writePositionsField( getInterfaceSize(), Foam::vector::zero );

        int globalOffset = 0;

        for ( int i = 0; i < Pstream::myProcNo(); i++ )
            globalOffset += nGlobalCenters[i];

        int offset = 0;

        forAll( movingPatchIDs, patchI )
        {
            const vectorField faceCentres( mesh.boundaryMesh()[movingPatchIDs[patchI]].faceCentres() );

            forAll( faceCentres, i )
            {
                writePositionsField[i + offset + globalOffset] = faceCentres[i];
            }

            offset += faceCentres.size();
        }

        reduce( writePositionsField, sumOp<vectorField>() );

        writePositions.resize( writePositionsField.size(), mesh.nGeometricD() );

        for ( int i = 0; i < writePositions.rows(); i++ )
            for ( int j = 0; j < writePositions.cols(); j++ )
                writePositions( i, j ) = writePositionsField[i][j];
    }
};

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

    assert( Pstream::nProcs() == 1 );

    std::shared_ptr<MapMeshSolidSolver> solid( new MapMeshSolidSolver( "solid", args, runTime ) );

    std::shared_ptr<foamFluidSolver> fluid( new CoupledFluidSolver( Foam::fvMesh::defaultRegion, args, runTime ) );

    /*
     * 1. Interpolate displacement from solid interface to fluid interface
     * 2. Mesh deformation fluid mesh
     */

    // Step 1: Interpolate displacement from solid interface to fluid interface

    std::shared_ptr<TPSFunction> tpsFunction( new TPSFunction() );
    rbf::RBFInterpolation rbfSolidToFluidInterface( tpsFunction );

    Eigen::MatrixXd positions, positionsInterpolation, values, valuesInterpolation;

    solid->getWritePositions( positions );
    fluid->getReadPositions( positionsInterpolation );

    rbfSolidToFluidInterface.compute( positions, positionsInterpolation );

    // valuesInterpolation = displacement of fluid interface
    // values = displacement of solid interface

    solid->getDisplacementLocal( values );

    rbfSolidToFluidInterface.interpolate( values, valuesInterpolation );

    // Step 2
    // Fluid mesh deformation

    fluid->setDisplacementLocal( valuesInterpolation );
    fluid->moveMesh();

    runTime->loop();
    runTime->write();
}
