
{
    phi = ( fvc::interpolate( HU ) / fvc::interpolate( AU ) ) & mesh.Sf();

    forAll( phi.boundaryField(), patchI )
    {
        if ( !phi.boundaryField()[patchI].coupled() )
        {
            phi.boundaryField()[patchI] =
                (
                U.boundaryField()[patchI]
                & mesh.Sf().boundaryField()[patchI]
                );
        }
    }

    phi += fvc::ddtPhiCorr( 1.0 / AU, U, phi );
}
