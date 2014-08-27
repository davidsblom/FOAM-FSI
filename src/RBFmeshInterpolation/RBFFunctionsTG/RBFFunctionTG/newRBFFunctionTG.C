
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFFunctionTG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<RBFFunctionTG> RBFFunctionTG::New
(
    const word& type,
    const dictionary& dict
)
{
    Info<< "Radial Basis Function interpolation: "
        << "Selecting RBF function: " << type << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "RBFFunctionTG::New(const word& type, const dictionary& dict)",
            dict
        )   << "Unknown RBFFunctionTG type "
            << type << endl << endl
            << "Valid  RBFFunctionTGs are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<RBFFunctionTG>(cstrIter()(dict.subDict(type + "Coeffs")));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
