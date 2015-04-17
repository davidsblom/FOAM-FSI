{
    Uf.oldTime();

    surfaceVectorField nf = mesh.Sf()/mesh.magSf();
    surfaceVectorField Utang = fvc::interpolate(U) - nf*(fvc::interpolate(U) & nf);
    surfaceVectorField Unor = phi/mesh.magSf() * nf;

    Uf = Utang + Unor;
}
