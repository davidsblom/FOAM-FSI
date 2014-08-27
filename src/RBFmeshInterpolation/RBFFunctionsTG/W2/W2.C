
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "W2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(W2, 0);
    addToRunTimeSelectionTable(RBFFunctionTG, W2, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::W2::W2(const scalar radius)
:
    RBFFunctionTG(),
    radius_(radius)
{}


Foam::W2::W2(const dictionary& dict)
:
    RBFFunctionTG(),
    radius_(readScalar(dict.lookup("radius")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::W2::~W2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::W2::weights
(
    const vectorField& points,
    const vector& controlPoint
) const
{
    scalarField dist = mag(points - controlPoint);

	//update by Thijs: Actually returning a tmp field
    tmp<scalarField> tRBF(new scalarField(dist.size(),0.0));
    scalarField& RBF = tRBF();
    //scalarField RBF(dist.size());

    RBF = neg(dist - radius_)*
        Foam::max(pow4(1 - (dist/radius_)), scalar(0))*(1 + 4*(dist/radius_));

    //update by Thijs: Actually returning a tmp field
	return tRBF;
    //return RBF;
}

Foam::tmp<Foam::scalarField> Foam::W2::weights
(
    const vectorField& points
) const
{
    tmp<scalarField> tRBF(new scalarField(points.size()*(points.size()+1)/2,0.0));
    scalarField& RBF = tRBF();

	for(int i = 0;i<points.size();i++)
	{
	    scalarField dist = mag(points - points[i]);
		for(int j = i;j<points.size();j++)
		{
		        RBF[i + (j+1)*j/2] = neg(dist[j] - radius_)*
        Foam::max(pow4(1 - (dist[j]/radius_)), scalar(0))*(1 + 4*(dist[j]/radius_));
		}
	}

	return tRBF;
}

// ************************************************************************* //
