
IOobject UfHeader
(
    "Uf",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

surfaceVectorField Uf
(
    UfHeader,
    mesh,
    dimensionedVector( "0", U.dimensions(), vector::zero )
);

if ( !UfHeader.headerOk() )
{
    Uf.oldTime();

    surfaceVectorField nf = mesh.Sf() / mesh.magSf();
    surfaceVectorField Utang = fvc::interpolate( U ) - nf * (fvc::interpolate( U ) & nf);
    surfaceVectorField Unor = phi / mesh.magSf() * nf;

    Uf = Utang + Unor;
}
