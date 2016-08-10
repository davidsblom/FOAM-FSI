
#include "fvCFD.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar ic_pressure_gauss ( scalar x, scalar y, scalar z)
{
    scalar pressure = 100000;
    //scalar d = ( x - 1 * x - 1 ) + y * y + z * z;
    scalar d = (x+1)*(x+1) + y * y + z * z;
    return pressure + 100 * std::exp ( -d / 0.2 * std::log( 2 ) );
}

int main(
    int argc,
    char * argv[]
    )
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Take the cell centres
    volVectorField centres = p.mesh().C();
    surfaceVectorField faceCentres = p.mesh().Cf();

    volScalarField x = centres.component(0);
    volScalarField y = centres.component(1);
    volScalarField z = centres.component(2);

    forAll( p.internalField(), i )
    {
        p.internalField()[i] = ic_pressure_gauss( x[i], y[i], z[i] );
    }

    forAll( p.boundaryField(), patchI )
    {
      const vectorField faceCentres( mesh.boundaryMesh()[patchI].faceCentres() );

      forAll( p.boundaryField()[patchI], i )
      {
          p.boundaryField()[patchI][i] = ic_pressure_gauss( faceCentres[i][0], faceCentres[i][1], faceCentres[i][2] );
      }
    }

    p.write();

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
