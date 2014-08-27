
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "Gauss.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Gauss, 0);
    addToRunTimeSelectionTable(RBFFunctionTG, Gauss, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Gauss::Gauss(const scalar radius)
:
    RBFFunctionTG(),
    radius_(radius)
{}


Foam::Gauss::Gauss(const dictionary& dict)
:
    RBFFunctionTG(),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Gauss::~Gauss()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::Gauss::weights
(
    const vectorField& points,
    const vector& controlPoint
) const
{
    // Algorithmic improvement, Matteo Lombardi.  21/Mar/2011
    scalarField sqrDist = magSqr(points - controlPoint);

//    return Foam::exp(-sqr(radius_)*sqrDist);
    return Foam::exp(-sqrDist/sqr(radius_));
}

Foam::tmp<Foam::scalarField> Foam::Gauss::weights
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
		        RBF[i + (j+1)*j/2] = Foam::exp(-sqrDist[j]/sqr(radius_));;
		}
	}

	return tRBF;
}

// ************************************************************************* //
