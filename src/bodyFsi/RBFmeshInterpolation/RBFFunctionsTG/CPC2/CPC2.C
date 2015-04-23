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

#include "CPC2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CPC2, 0);
    addToRunTimeSelectionTable(RBFFunctionTG, CPC2, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CPC2::CPC2(const scalar radius)
:
    RBFFunctionTG(),
    radius_(radius)
{}


Foam::CPC2::CPC2(const dictionary& dict)
:
    RBFFunctionTG(),
    radius_(readScalar(dict.lookup("radius")))
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CPC2::~CPC2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::CPC2::weights
(
    const vectorField& points,
    const vector& controlPoint
) const
{
    scalarField dist = mag(points - controlPoint)/radius_;

	//update by Thijs: Actually returning a tmp field
    tmp<scalarField> tRBF(new scalarField(dist.size(),0.0));
    scalarField& RBF = tRBF();

    forAll(RBF, i)
    {
        if (dist[i] < 1.0)
        {
            RBF[i] = (1-dist[i])*(1-dist[i])*(1-dist[i])*(1-dist[i]) * (4*dist[i]+1);
        }
        else
        {
            RBF[i] = 0.0;
        }
    }
    
	//update by Thijs: Actually returning a tmp field
	return tRBF;
    //return RBF;
}

Foam::tmp<Foam::scalarField> Foam::CPC2::weights
(
    const vectorField& points
) const
{
    tmp<scalarField> tRBF(new scalarField(points.size()*(points.size()+1)/2,0.0));
    scalarField& RBF = tRBF();

	for(int i = 0;i<points.size();i++)
	{
	    scalarField dist = mag(points - points[i])/radius_;
		for(int j = i;j<points.size();j++)
		{
			if (dist[j] < 1.0)
		    {
                RBF[i + (j+1)*j/2] = (1-dist[j])*(1-dist[j])*(1-dist[j])*(1-dist[j]) * (4*dist[j]+1);
		    }
		    else
		    {
		        RBF[i + (j+1)*j/2] = 0.0;
		    }
		}
	}

	return tRBF;
}


// ************************************************************************* //
