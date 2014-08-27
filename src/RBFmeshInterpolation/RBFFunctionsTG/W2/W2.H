
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef W2_H
#define W2_H

#include "RBFFunctionTG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              Class W2 Declaration
\*---------------------------------------------------------------------------*/

class W2
:
    public RBFFunctionTG
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        W2(const W2&);

        //- Disallow default bitwise assignment
        void operator=(const W2&);


public:

    //- Runtime type information
    TypeName("W2");

    // Constructors

        //- Construct given radius
        W2(const scalar radius);

        //- Construct from dictionary
        W2(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunctionTG> clone() const
        {
            return autoPtr<RBFFunctionTG>(new W2(this->radius_));
        }


    // Destructor

        virtual ~W2();


    // Member Functions

        //- Return weights given points
        virtual tmp<scalarField> weights
        (
            const vectorField& points,
            const vector& controlPoint
        ) const;

        virtual tmp<scalarField> weights
        (
        	const vectorField& points
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
