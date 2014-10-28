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

Class
    solidPolyMesh

Description
    Mesh for solid mechanics allowing implicit cell point neighbours
    to enable a block coupled solution.

    Note on addressing:
        The problem of addressing with fvMesh is that only cell face neighbours
        are treated implicitly; however, for a fully coupled solution in solid
        mechanics all cell point neighbours must be implicit. Therefore, we must
        create a new lduAddressing to allow point neighbours.

        To achieve this we will define imaginary edges which connect cell
        centres to the cell centres of point neighbours. This way the imaginary
        edges act in a analogous manner to the faces in the fvMesh, where each
        imaginary edge has an owner and neighbour. This allows us to create
        lduAddressing in a similar manner to fvMesh; the size of upper and lower
        will be numImagEdges. The imagEdges are created by visiting every cell
        in ascending order and adding all point neighbour in ascending order not
        already added. This way the starts of the imagEdges define the owners
        (and lowerAddr) and the ends define the neighbours (and upperAddr).

SourceFiles
    solidPolyMesh.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef solidPolyMesh_H
#define solidPolyMesh_H

#include "GeoMesh.H"
#include "polyMesh.H"
#include "lduMesh.H"
#include "primitiveMesh.H"
#include "solidPolyBoundaryMesh.H"
#include "solidPolyMeshLduAddressing.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solidPolyMeshLduAddressing;

/*---------------------------------------------------------------------------*\
                     Class solidPolyMesh Declaration
\*---------------------------------------------------------------------------*/

class solidPolyMesh
:
    public GeoMesh<polyMesh>,
    public lduMesh
{
    // Permanent data

        //- Boundary mesh
        solidPolyBoundaryMesh boundary_;


    // Demand-driven data

        //- LDU addressing
        mutable solidPolyMeshLduAddressing* lduPtr_;

        //- Imaginary edges
        mutable edgeList* imagEdgesPtr_;

        //- Whether the imaginary edge is near (ie connected to a near
        // neighbour) or not
        mutable boolList* nearImagEdgesPtr_;

        //- This map tranforms the fvMesh upper/lower to the new addressing
        mutable labelList* fvMeshAddressingMapPtr_;

        //- Interpolation scheme weighting factor array for interpolation from
        // cell centres to points
        mutable scalarListList* pointWeightsPtr_;

        //- Interpolation scheme weighting factor array for interpolation from
        // cell centres to edge centres
        mutable scalarListList* edgeWeightsPtr_;

        //- Cell point cells i.e. cell neighbours connect through any of the
        // cell points
        mutable labelListList* cellPointCellsPtr_;

        //- Cell imaginary edges
        mutable labelListList* cellImagEdgesPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        solidPolyMesh(const solidPolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const solidPolyMesh&);


    // Private member functions to calculate demand driven data

        //- Calculate imaginary edges and set nearImagEdges
        void calcImagEdges() const;

        //- Calculate fvMesh addressing map
        void calcFvMeshAddressingMap() const;

        //- Calculte interpolation scheme weighting factor array
        void makePointWeights() const;

        //- Calculte interpolation scheme weighting factor array
        void makeEdgeWeights() const;

        //- Calculte cell imaginary edges
        void calcCellImagEdges() const;

        //- Clear out all data
        void clearOut() const;


public:

    // Declare name of the class and its debug switch
    ClassName("solidPolyMesh");

    typedef solidPolyMesh Mesh;
    typedef solidPolyBoundaryMesh BoundaryMesh;

        // Static mesh analysis data

             //- Estimated number of imaginary edges per cell
             static const unsigned imagEdgesPerCell_ = 40;


    // Constructors

        //- Construct from components
        explicit solidPolyMesh(const polyMesh& pMesh);

    // Destructor

        ~solidPolyMesh();


    // Member Functions

        // Access

            //- Return the top-level database
            const Time& time() const
            {
                return mesh_.time();
            }

            //- Return the object registry
            virtual const objectRegistry& thisDb() const
            {
                return GeoMesh<polyMesh>::thisDb();
            }

            //- Return imaginary edges
            const edgeList& imagEdges() const;

            //- Return near imaginary edges
            const boolList& nearImagEdges() const;

            //- Return fvMesh addressing map
            const labelList& fvMeshAddressingMap() const;

            //- Return interpolation point cell weighting array
            const scalarListList& pointWeights() const;

            //- Return interpolation edge cell weighting array
            const scalarListList& edgeWeights() const;

            //- Return cell point cells
            const labelListList& cellPointCells() const;

            //- Return cell imaginary edges
            const labelListList& cellImagEdges() const;

            //- Return number of cells in polyhedral mesh
            label nCells() const
            {
                return mesh_.nCells();
            }

            //- Return reference to boundary mesh
            const solidPolyBoundaryMesh& boundary() const
            {
                 return boundary_;
            }

            //- Return ldu addressing
            const lduAddressing& lduAddr() const;

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            virtual lduInterfacePtrsList interfaces() const
            {
                return boundary().interfaces();
            }

            //- Owner
            const unallocLabelList& owner() const
            {
                return lduAddr().lowerAddr();
            }

            //- Neighbour
            const unallocLabelList& neighbour() const
            {
                return lduAddr().upperAddr();
            }

            //- Return imaginary edge ID given two cell IDs
            // Retruns -1 if the two cells are not connected by an imaginary
            // edge
            label findImaginaryEdge
            (
                const label firstCellID, const label secondCellID
            ) const;

            // Parallel mesh analysis data

                //- Return parallel info
                const globalMeshData& globalData() const
                {
                    return mesh_.globalData();
                }

                //- Shared parallel points
                //const labelList& parallelPoints() const;

                //- Shared parallel edges
                // const edgeList& parallelEdges() const;


        // Edit

            //- Update mesh topology using the morph engine
            // void updateMesh(const solidPolyMeshMapper& mapper);


    // Member Operators

        bool operator!=(const solidPolyMesh&) const;
        bool operator==(const solidPolyMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
