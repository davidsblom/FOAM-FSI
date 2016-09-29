
#include "fvCFD.H"
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar ic_pressure_gauss ( scalar x, scalar y, scalar z)
{
    scalar pressure = 100000;
    //scalar d = ( x - 1 * x - 1 ) + y * y + z * z;
    scalar d = (x+1)*(x+1) + y * y + z * z;
    return pressure + 50 * std::exp ( -d / 0.1 * std::log( 2 ) );
}

scalar ic_temperature_gauss ( scalar x, scalar y, scalar z)
{
    scalar temperature = 275.78599007170436;
    //scalar d = ( x - 1 * x - 1 ) + y * y + z * z;
    scalar d = (x+1)*(x+1) + y * y + z * z;
    return temperature + 25 * std::exp ( -d / 0.1 * std::log( 2 ) );
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

    // volScalarField p
    // (
    //     IOobject
    //     (
    //         "p",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::NO_WRITE
    //     ),
    //     mesh
    // );

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // // Take the cell centres
    // volVectorField centres = p.mesh().C();
    // surfaceVectorField faceCentres = p.mesh().Cf();
    //
    // volScalarField x = centres.component(0);
    // volScalarField y = centres.component(1);
    // volScalarField z = centres.component(2);
    //
    // forAll( p.internalField(), i )
    // {
    //     p.internalField()[i] = ic_pressure_gauss( x[i], y[i], z[i] );
    // }
    //
    // forAll( p.boundaryField(), patchI )
    // {
    //   const vectorField faceCentres( mesh.boundaryMesh()[patchI].faceCentres() );
    //
    //   forAll( p.boundaryField()[patchI], i )
    //   {
    //       p.boundaryField()[patchI][i] = ic_pressure_gauss( faceCentres[i][0], faceCentres[i][1], faceCentres[i][2] );
    //   }
    // }
    //
    // p.write();

    // Take the cell centres
    volVectorField centres = T.mesh().C();
    surfaceVectorField faceCentres = T.mesh().Cf();

    volScalarField x = centres.component(0);
    volScalarField y = centres.component(1);
    volScalarField z = centres.component(2);

    forAll( T.internalField(), i )
    {
        T.internalField()[i] = ic_temperature_gauss( x[i], y[i], z[i] );
    }

    forAll( T.boundaryField(), patchI )
    {
      const vectorField faceCentres( mesh.boundaryMesh()[patchI].faceCentres() );

      forAll( T.boundaryField()[patchI], i )
      {
          T.boundaryField()[patchI][i] = ic_temperature_gauss( faceCentres[i][0], faceCentres[i][1], faceCentres[i][2] );
      }
    }

    T.write();

    Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << nl << endl;

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
