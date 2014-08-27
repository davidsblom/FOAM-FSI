
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFFunctionTG_H
#define RBFFunctionTG_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "tmp.H"
#include "autoPtr.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RBFFunctionTG Declaration
\*---------------------------------------------------------------------------*/

class RBFFunctionTG
{
    // Private Member Functions

        //- Disallow copy construct
        RBFFunctionTG(const RBFFunctionTG&);

        //- Disallow default bitwise assignment
        void operator=(const RBFFunctionTG&);


public:

    //- Runtime type information
    TypeName("RBFFunctionTG");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RBFFunctionTG,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );


    // Selectors

        //- Return a pointer to the selected RBF function
        static autoPtr<RBFFunctionTG> New
        (
            const word& type,
            const dictionary& dict
        );


    // Constructors

        //- Construct null
        RBFFunctionTG()
        {}

        //- Create and return a clone
        virtual autoPtr<RBFFunctionTG> clone() const = 0;


    // Destructor

        virtual ~RBFFunctionTG()
        {}


    // Member Functions

        //- Return RBF weights
        virtual tmp<scalarField> weights
        (
            const vectorField& points,
            const vector& controlPoint
        ) const = 0;

        //- Return RBF weights for control matrix
        virtual tmp<scalarField> weights
        (
        	const vectorField& points
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
