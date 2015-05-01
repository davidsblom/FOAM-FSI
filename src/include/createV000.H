
// Create object to hold V000
IOobject V000Header
(
    "V000",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

volScalarField V000
(
    V000Header,
    mesh,
    dimensionedScalar( "0", dimVolume, 0.0 )
);

// Check if scheme is bdf3. If not do not fill V000
const IOdictionary & fvSchemes = mesh.lookupObject<IOdictionary>( "fvSchemes" );
const dictionary & ddtSchemes = fvSchemes.subDict( "ddtSchemes" );

word ddtScheme;

if ( ddtSchemes.found( "ddt(U)" ) )
{
    ddtScheme = word( ddtSchemes.lookup( "ddt(U)" ) );
}
else
{
    ddtScheme = word( ddtSchemes.lookup( "default" ) );
}

bool updateV000 = false;
bool firstV000 = true;

if ( ddtScheme == word( "bdf3" ) )
{
    updateV000 = true;

    if ( !V000Header.headerOk() )
    {
        V000.internalField() = mesh.V();
    }
}
