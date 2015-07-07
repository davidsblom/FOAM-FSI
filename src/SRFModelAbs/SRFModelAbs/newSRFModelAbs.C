
#include "SRFModelAbs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace SRF
    {
        // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

        autoPtr<SRFModelAbs> SRFModelAbs::New
            ( const volVectorField & U )
        {
            word SRFModelAbsTypeName;

            // Enclose the creation of the SRFPropertiesDict to ensure it is
            // deleted before the SRFModelAbs is created - otherwise the dictionary
            // is entered in the database twice
            {
                IOdictionary SRFPropertiesDict
                (
                    IOobject
                    (
                        "SRFProperties",
                        U.time().constant(),
                        U.db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );

                SRFPropertiesDict.lookup( "SRFModelAbs" ) >> SRFModelAbsTypeName;
            }

            Info << "Selecting SRFModelAbs " << SRFModelAbsTypeName << endl;

            dictionaryConstructorTable::iterator cstrIter =
                dictionaryConstructorTablePtr_->find( SRFModelAbsTypeName );

            if ( cstrIter == dictionaryConstructorTablePtr_->end() )
            {
                FatalErrorIn
                (
                    "SRFModelAbs::New(const fvMesh&)"
                ) << "Unknown SRFModelAbs type " << SRFModelAbsTypeName
                  << nl << nl
                  << "Valid SRFModelAbs types are :" << nl
                  << dictionaryConstructorTablePtr_->sortedToc()
                  << exit( FatalError );
            }

            return autoPtr<SRFModelAbs>( cstrIter() ( U ) );
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    } // End namespace SRF
} // End namespace Foam

// ************************************************************************* //
