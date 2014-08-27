
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Gauss_H
#define Gauss_H

#include "RBFFunctionTG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              Class Gauss Declaration
\*---------------------------------------------------------------------------*/

class Gauss
:
    public RBFFunctionTG
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        Gauss(const Gauss&);

        //- Disallow default bitwise assignment
        void operator=(const Gauss&);


public:

    //- Runtime type information
    TypeName("Gauss");

    // Constructors

        //- Construct given radius
        Gauss(const scalar radius);

        //- Construct from dictionary
        Gauss(const dictionary& dict);

        //- Create and return a clone
        virtual autoPtr<RBFFunctionTG> clone() const
        {
            return autoPtr<RBFFunctionTG>(new Gauss(this->radius_));
        }


    // Destructor

        virtual ~Gauss();


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
