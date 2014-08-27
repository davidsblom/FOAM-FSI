
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef CPC2_H
#define CPC2_H

#include "RBFFunctionTG.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                              Class CPC2 Declaration
\*---------------------------------------------------------------------------*/

class CPC2
:
    public RBFFunctionTG
{
    // Private data

        //- Radius
        scalar radius_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CPC2(const CPC2&);

        //- Disallow default bitwise assignment
        void operator=(const CPC2&);


public:

    //- Runtime type information
    TypeName("CPC2");

    // Constructors

        //- Construct given radius
        CPC2(const scalar radius);

        //- Construct from dictionary
        CPC2(const dictionary& dict);

        virtual autoPtr<RBFFunctionTG> clone() const
        {
            return autoPtr<RBFFunctionTG>(new CPC2(this->radius_));
        }


    // Destructor

        virtual ~CPC2();


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
