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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::aeroForces

Description
    Calculates the forces and moments by integrating the pressure and
    skin-friction forces over a given list of patches.

    Member function calcForcesMoment()calculates and returns the forces and
    moments.

    Member function forces::write() calls calcForcesMoment() and writes the
    forces and moments into the file \<timeDir\>/forces.dat

SourceFiles
    forces.C

\*---------------------------------------------------------------------------*/

#ifndef aeroForces_H
#define aeroForces_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class forces Declaration
\*---------------------------------------------------------------------------*/

class aeroForces
{
protected:

    // Private data

        //- Name of this set of forces,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        // Read from dictionary

            //- Patches to integrate forces over
            const labelList& patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Name of density field (optional)
            word rhoName_;

            //- Reference density needed for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;

        //- Forces/moment file ptr
        autoPtr<OFstream> forcesFilePtr_;

    // Private Member Functions

        //- If the forces file has not been created create it
        void makeFile();

        //- Output file header information
        void writeFileHeader();

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devRhoReff() const;

        //- Return rho if rhoName is specified otherwise rhoRef
        tmp<volScalarField> rho() const;

        //- Return rhoRef if the pressure field is dynamic, i.e. p/rho
        //  otherwise return 1
        scalar rho(const volScalarField& p) const;

        //- Disallow default bitwise copy construct
        aeroForces(const aeroForces&);

        //- Disallow default bitwise assignment
        void operator=(const aeroForces&);


public:

    //- Runtime type information
    TypeName("aeroForces");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        aeroForces
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const labelList& movingPatchIDs
        );


    //- Destructor
    virtual ~aeroForces();

    // Member Functions

        //- Return name of the set of forces
        const word& name() const
        {
            return name_;
        }

        //- Read the forces data
        void read(const dictionary&);

        //- Write the forces
        void write();

        Field<vectorField> calcForceFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //