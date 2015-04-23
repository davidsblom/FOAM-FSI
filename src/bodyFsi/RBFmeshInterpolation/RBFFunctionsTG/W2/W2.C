/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

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
