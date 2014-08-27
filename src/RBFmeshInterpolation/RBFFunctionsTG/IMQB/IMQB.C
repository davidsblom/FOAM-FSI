
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "IMQB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IMQB, 0);
    addToRunTimeSelectionTable(RBFFunctionTG, IMQB, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IMQB::IMQB(const scalar radius)
:
    RBFFunctionTG(),
    radius_(radius)
{}


// Construct from dictionary
Foam::IMQB::IMQB(const dictionary& dict)
:
    RBFFunctionTG(),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::IMQB::~IMQB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::IMQB::weights
(
    const vectorField& points,
    const vector& controlPoint
) const
{
    // Algorithmic improvement, Matteo Lombardi.  21/Mar/2011
    scalarField sqrDist = magSqr(points - controlPoint);

    return 1/sqrt(sqrDist + sqr(radius_));
}

Foam::tmp<Foam::scalarField> Foam::IMQB::weights
(
    const vectorField& points
) const
{
    tmp<scalarField> tRBF(new scalarField(points.size()*(points.size()+1)/2,0.0));
    scalarField& RBF = tRBF();

	for(int i = 0;i<points.size();i++)
	{
	    scalarField sqrDist = magSqr(points - points[i]);
		for(int j = i;j<points.size();j++)
		{
		        RBF[i + (j+1)*j/2] = 1/sqrt(sqrDist[j] + sqr(radius_));
		}
	}

	return tRBF;
}


// ************************************************************************* //
