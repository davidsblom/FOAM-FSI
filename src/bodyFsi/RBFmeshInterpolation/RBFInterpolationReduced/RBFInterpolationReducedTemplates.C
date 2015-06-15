/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    RBF interpolation templates

Author
    Thijs Gillebaart, TU Delft

\*---------------------------------------------------------------------------*/

#include "RBFInterpolationReduced.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::RBFInterpolationReduced::interpolateRed
(
  const Field<Type>& ctrlField
) const  
{
    if (ctrlField.size() != NRC)
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > RBFInterpolationReduced::interpolateRed\n"
            "(\n"
            "    const Field<Type>& ctrlField\n"
            ") const"
        )   << "Incorrect size of source field.  Size = " << ctrlField.size()
            << " nControlPoints = " << NRC
            << abort(FatalError);
    }

    //set dimensions
    int dimtmp = 3;
    if(Dimension_ == "TwoD"){
        Info << "Just doing 2D RBF calulations" << endl;
        dimtmp = 2;
    }else if(Dimension_ == "ThreeD"){
        Info << "Doing full blown 3D RBF calculations" << endl;
        dimtmp = 3;
    }else{
        FatalErrorIn("void Foam::RBFInterpolationReduced::createReducedEvaluationMatrix(labelList& movingControlIndex)")
            << "Dimension are neither TwoD or ThreeD"
            << abort(FatalError);
    }
    const int DIM = dimtmp;

    tmp<Field<Type> > tresult
    (
        new Field<Type>(internalPoints_.size(), pTraits<Type>::zero) 
    );
    Field<Type>& result = tresult(); 

    // Construct dmc
    double* dmc = new double[NRC*DIM];
    int i,j,k;
    for ( i = 0; i < NRC; i++ )    
    {
        for( k = 0; k < DIM ; k++ )
        {
            dmc[i+k*NRC] = ctrlField[i][k];
        }
    }

    //Calculate internal field displacement for each point
    for ( i =0; i < internalPoints_.size(); i++ )
    {
        double* Hi = new double[NRC];

        for ( j=0; j<NRC; j++ )
        {
            Hi[j]     = Hred[i+j*internalPoints_.size()];  
        }

        // initializing output 
        double* Result = new double[DIM];
        for ( j = 0; j < DIM; j++ ) 
        {
            Result[j] = 0.0;
        }

        // coefficients for solving matrix-vector product with BLAS2
        const double alpha1  = 1;
        const double alpha2  = 0;
        int incx = 1;
        int incy = 1;

        //NRHS: http://for2html.sourceforge.net/BLAS/dgemv.html
        dgemv_ ( "T", &NRC, &DIM, &alpha1, dmc, &NRC, Hi, &incx, &alpha2, Result, &incy );

        for (k = 0; k < DIM ; k++)
        {
            result[i][k] = Result[k];
        }

        //clean up
        delete[] Hi;
        delete[] Result;
    }
    delete[] dmc;

    //Return result
    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //