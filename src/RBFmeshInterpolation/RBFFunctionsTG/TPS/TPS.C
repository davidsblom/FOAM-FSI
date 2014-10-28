
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "TPS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug( TPS, 0 );
  addToRunTimeSelectionTable( RBFFunctionTG, TPS, dictionary );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TPS::TPS( const scalar radius )
  :
  RBFFunctionTG(),
  radius_( radius )
{}


Foam::TPS::TPS( const dictionary & dict )
  :
  RBFFunctionTG(),
  radius_( readScalar( dict.lookup( "radius" ) ) )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TPS::~TPS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::TPS::weights
(
  const vectorField & points,
  const vector & controlPoint
) const
{
  scalarField dist = mag( points - controlPoint );

  // update by Thijs: Actually returning a tmp field
  tmp<scalarField> tRBF( new scalarField( dist.size(), 0.0 ) );
  scalarField & RBF = tRBF();

  forAll( RBF, i )
  {
    if ( dist[i] > SMALL )
    {
      RBF[i] = sqr( dist[i] ) * log( dist[i] );
    }
    else
    {
      RBF[i] = 0.0;
    }
  }

  /*forAll(RBF, i)
   * {
   *  RBF[i] = sqr(dist[i])*log(dist[i]);
   * }*/

  // update by Thijs: Actually returning a tmp field
  return tRBF;
}

Foam::tmp<Foam::scalarField> Foam::TPS::weights
  ( const vectorField & points ) const
{
  tmp<scalarField> tRBF( new scalarField( points.size() * (points.size() + 1) / 2, 0.0 ) );
  scalarField & RBF = tRBF();

  for ( int i = 0; i < points.size(); i++ )
  {
    scalarField dist = mag( points - points[i] );

    for ( int j = i; j < points.size(); j++ )
    {
      if ( dist[j] > SMALL )
      {
        RBF[i + (j + 1) * j / 2] = sqr( dist[j] ) * log( dist[j] );
      }
      else
      {
        RBF[i + (j + 1) * j / 2] = 0.0;
      }
    }
  }

  return tRBF;
}

// ************************************************************************* //
