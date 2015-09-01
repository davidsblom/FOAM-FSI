
#   include "readTimeControls.H"
#   include "readAutoPIMPLEControls.H"

bool correctPhi = false;

if ( pimple.found( "correctPhi" ) )
{
    correctPhi = Switch( pimple.lookup( "correctPhi" ) );
}

bool checkMeshCourantNo = false;

if ( pimple.found( "checkMeshCourantNo" ) )
{
    checkMeshCourantNo = Switch( pimple.lookup( "checkMeshCourantNo" ) );
}
