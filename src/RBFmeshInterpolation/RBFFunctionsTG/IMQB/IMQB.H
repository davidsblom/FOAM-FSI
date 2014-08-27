
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef IMQB_H
#define IMQB_H

#include "RBFFunctionTG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              Class IMQB Declaration
\*---------------------------------------------------------------------------*/

class IMQB
:
    public RBFFunctionTG
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IMQB(const IMQB&);

        //- Disallow default bitwise assignment
        void operator=(const IMQB&);


public:

    //- Runtime type information
    TypeName("IMQB");

    // Constructors

        //- Construct given radius
        IMQB(const scalar radius);

        //- Construct from dictionary
        IMQB(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunctionTG> clone() const
        {
            return autoPtr<RBFFunctionTG>(new IMQB(this->radius_));
        }


    // Destructor

        virtual ~IMQB();


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
