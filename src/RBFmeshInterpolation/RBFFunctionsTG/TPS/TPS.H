
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TPS_H
#define TPS_H

#include "RBFFunctionTG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              Class TPS Declaration
\*---------------------------------------------------------------------------*/

class TPS
:
    public RBFFunctionTG
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        TPS(const TPS&);

        //- Disallow default bitwise assignment
        void operator=(const TPS&);


public:

    //- Runtime type information
    TypeName("TPS");

    // Constructors

        //- Construct given radius
        TPS(const scalar radius);

        //- Construct from dictionary
        TPS(const dictionary& dict);

        virtual autoPtr<RBFFunctionTG> clone() const
        {
            return autoPtr<RBFFunctionTG>(new TPS(this->radius_));
        }


    // Destructor

        virtual ~TPS();


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
