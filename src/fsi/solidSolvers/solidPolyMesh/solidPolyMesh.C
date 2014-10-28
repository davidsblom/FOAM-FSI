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
    The tetrahedral decomposition of the mesh based on the underlying polyMesh.
    Uses minimal storage - the tet decomposition is provided on a cell-by-cell
    basis

\*---------------------------------------------------------------------------*/

#include "solidPolyMesh.H"
#include "demandDrivenData.H"
#include "SortableList.H"

#include "pointFields.H"
#include "coupledPointPatchFields.H"
#include "pointConstraint.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::solidPolyMesh, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void solidPolyMesh::clearOut() const
{
    deleteDemandDrivenData(lduPtr_);
    deleteDemandDrivenData(imagEdgesPtr_);
    deleteDemandDrivenData(nearImagEdgesPtr_);
    deleteDemandDrivenData(fvMeshAddressingMapPtr_);
    deleteDemandDrivenData(pointWeightsPtr_);
    deleteDemandDrivenData(edgeWeightsPtr_);
    deleteDemandDrivenData(cellPointCellsPtr_);
    deleteDemandDrivenData(cellImagEdgesPtr_);
}


void solidPolyMesh::calcImagEdges() const
{
    if (debug)
    {
        Info<< "void solidPolyMesh::calcImagEdges() const "
            << "creating imaginary edges" << endl;
    }

    // It is an error to attempt to recalculate if the pointer is already set
    if (imagEdgesPtr_ || nearImagEdgesPtr_)
    {
        FatalErrorIn("solidPolyMesh::calcImagEdges() const")
            << "imaginary edges, near imaginary edges, or cell point cells"
            << "lists already exist" << abort(FatalError);
    }

    // Define imaginery edges connecting cell centres to the cell centres of
    // surrounding point-cell neighbours
    // We will follow a similar procedure to primitiveMesh::calcEdges() to
    // calculate the edges.

    // Count the maximun number of imaginary edges
    label maxEdges = 0;

    const polyMesh& mesh = (*this)();

    const cellList& c = mesh.cells();
    const labelListList& pc = mesh.pointCells();
    const labelListList& cp = mesh.cellPoints();
    const labelListList& cc = mesh.cellCells();

    // We will create the list of imaginary neighbours
    // Imaginary neighbours are the neighbours cell connected through a point
    labelListList imagNeiOfCell(c.size());

    // Create cellPointCells array
    cellPointCellsPtr_ = new labelListList(mesh.nCells());
    labelListList& cellPointCells = *cellPointCellsPtr_;

    // Estimate of number of imagEdges per cell
    //const label edgesPerCell = 26;

    forAll(c, cellI)
    {
        // Count cell imagEdges i.e. cell point cells and remember them
        label ie = 0;
        labelHashSet curCellPointCells(imagEdgesPerCell_);
        forAll(cp[cellI], cellPointI)
        {
            // Count point cells for every point in the cell. The minus 1 is
            // because we don't count the current cell i.e. if there are 4
            // pointCells then we have 3 imaginary edges from our original
            // cell
            forAll(pc[cp[cellI][cellPointI]], pointCellI)
            {
                // Add the cell if it is not the original cell and it has
                // not already been added
                const label otherCellID =
                    pc[cp[cellI][cellPointI]][pointCellI];
                if (otherCellID != cellI)
                {
                    if (!curCellPointCells.found(otherCellID))
                    {
                        curCellPointCells.insert(otherCellID);
                        ie++;
                    }
                }
            }
        }

        maxEdges += ie;

        // Record imaginary neighbours using the cellPointCells
        edgeList cellImagEdges(ie);
        labelList cellImagNeiCells(ie);
        //const labelList cpc = curCellPointCells.toc();
        labelList& cpc = cellPointCells[cellI];
        cpc.setSize(curCellPointCells.size());
        cpc = curCellPointCells.toc();

        forAll(cellImagEdges, edgeI)
        {
            label otherCellID = cpc[edgeI];
            if (otherCellID != cellI)
            {
                cellImagNeiCells[edgeI] = otherCellID;
            }
        }

        imagNeiOfCell[cellI] = cellImagNeiCells;
    }

    // Edge calculation
    // The imagEdges will be ordered in an analagous way to the faces in a
    // standard fvMesh i.e. we start at cell 0 and add all imagEdges, then
    // go to cell 1 and add all imagEdges not already added, the cell 2 etc.
    // Within a cell we add the imagEdges in order of increasing neighbour
    // cell ID
    imagEdgesPtr_ = new edgeList(maxEdges);
    edgeList& imagEdges = *imagEdgesPtr_;
    nearImagEdgesPtr_ = new boolList(maxEdges, false);
    boolList& nearImagEdges = *nearImagEdgesPtr_;

    label nEdges = 0;

    forAll(c, cellI)
    {
        const labelList cellCells = cc[cellI];

        // Sort the imagNei of this cell
        SortableList<label> sortedImagNei(imagNeiOfCell[cellI]);
        forAll(sortedImagNei, imagNeiCellI)
        {
            const label otherCellID = sortedImagNei[imagNeiCellI];
            if (otherCellID > cellI)
            {
                imagEdges[nEdges] = edge(cellI, otherCellID);

                // Check if the other cell is a near neighbour
                forAll(cellCells, cellCellI)
                {
                    if (otherCellID == cellCells[cellCellI])
                    {
                        nearImagEdges[nEdges] = true;
                        break;
                    }
                }
                nEdges++;
            }
        }
    }

    // Resize to actual number of imagEdges
    imagEdges.resize(nEdges);
    nearImagEdges.resize(nEdges);
}


void solidPolyMesh::calcFvMeshAddressingMap() const
{
    const polyMesh& mesh = (*this)();

    fvMeshAddressingMapPtr_ = new labelList(mesh.nInternalFaces(), -1);
    labelList& fvMap = *fvMeshAddressingMapPtr_;

    const boolList& nearImEdges = nearImagEdges();

    label fvFaceI = 0;
    forAll(nearImEdges, imEdgeI)
    {
        if (nearImEdges[imEdgeI])
        {
            fvMap[fvFaceI] = imEdgeI;
            fvFaceI++;
        }
    }

    if (min(fvMap) < 0)
    {
        FatalErrorIn
        (
            "void solidPolyMesh::calcFvMeshAddressingMap() const"
        )   << "All faces not set" << abort(FatalError);
    }
}


void solidPolyMesh::makePointWeights() const
{
    // Taken from volPointInterpolation.C

    if (debug)
    {
        Info<< "solidPolyMesh::makePointWeights() : "
            << "constructing weighting factors"
            << endl;
    }

    const polyMesh& mesh = (*this)();

    const pointField& points = mesh.points();
    const labelListList& pointCells = mesh.pointCells();
    const vectorField& cellCentres = mesh.cellCentres();

    // Allocate storage for weighting factors
    pointWeightsPtr_ = new scalarListList(points.size());
    scalarListList& pointWeights = *pointWeightsPtr_;
    //pointWeights_.clear();
    //pointWeights_.setSize(points.size());

    forAll(pointWeights, pointi)
    {
        pointWeights[pointi].setSize(pointCells[pointi].size());
    }

    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh.polyMesh::instance(),
            mesh
        ),
        pointMesh::New(mesh),
        dimensionedScalar("zero", dimless, 0)
    );

    // Calculate inverse distances between cell centres and points
    // and store in weighting factor array
    forAll(points, pointi)
    {
        scalarList& pw = pointWeights[pointi];
        const labelList& pcp = pointCells[pointi];

        forAll(pcp, pointCelli)
        {
            pw[pointCelli] =
                1/mag(points[pointi] - cellCentres[pcp[pointCelli]]);

            sumWeights[pointi] += pw[pointCelli];
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].initAddField();
        }
    }

    forAll(sumWeights.boundaryField(), patchi)
    {
        if (sumWeights.boundaryField()[patchi].coupled())
        {
            sumWeights.boundaryField()[patchi].addField
            (
                sumWeights.internalField()
            );
        }
    }

    forAll(points, pointi)
    {
        scalarList& pw = pointWeights[pointi];

        forAll(pw, pointCelli)
        {
            pw[pointCelli] /= sumWeights[pointi];
        }
    }

    if (debug)
    {
        Info<< "solidPolyMesh::makePointWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


void solidPolyMesh::makeEdgeWeights() const
{
    if (debug)
    {
        Info<< "solidPolyMesh::makeEdgeWeights() : "
            << "constructing edge weighting factors"
            << endl;
    }

    const polyMesh& mesh = (*this)();

    const pointField& points = mesh.points();
    const edgeList edges = mesh.edges();
    const labelListList& edgeCells = mesh.edgeCells();
    const vectorField& cellCentres = mesh.cellCentres();

    // Allocate storage for weighting factors
    edgeWeightsPtr_ = new scalarListList(edges.size());
    scalarListList& edgeWeights = *edgeWeightsPtr_;

    forAll(edgeWeights, edgei)
    {
        edgeWeights[edgei].setSize(edgeCells[edgei].size());
    }

    scalarField sumWeights(edgeWeights.size(), 0.0);

    // Calculate inverse distances between cell centres and edge centres
    // and store in weighting factor array
    forAll(edges, edgei)
    {
        scalarList& ew = edgeWeights[edgei];
        const labelList& ec = edgeCells[edgei];

        forAll(ec, edgeCelli)
        {
            vector edgeCentre = edges[edgei].centre(points);

            ew[edgeCelli] =
                1/mag(edgeCentre - cellCentres[ec[edgeCelli]]);

            sumWeights[edgei] += ew[edgeCelli];
        }
    }

    if (Pstream::parRun())
    {
        WarningIn
        (
            "solidPolyMesh::makeEdgeWeights"
        ) << "The weights are currently not becorrect on processor boundaries"
          //<< abort(FatalError);
          << endl;
    }

    // Parallel method for point fields
    // Similar method can be used for edge field manually
    // forAll(sumWeights.boundaryField(), patchi)
    // {
    //     if (sumWeights.boundaryField()[patchi].coupled())
    //     {
    //         sumWeights.boundaryField()[patchi].initAddField();
    //     }
    // }

    // forAll(sumWeights.boundaryField(), patchi)
    // {
    //     if (sumWeights.boundaryField()[patchi].coupled())
    //     {
    //         sumWeights.boundaryField()[patchi].addField
    //         (
    //             sumWeights.internalField()
    //         );
    //     }
    // }

    forAll(edges, edgei)
    {
        scalarList& pw = edgeWeights[edgei];

        forAll(pw, edgeCelli)
        {
            pw[edgeCelli] /= sumWeights[edgei];
        }
    }

    if (debug)
    {
        Info<< "solidPolyMesh::makeEdgeWeights() : "
            << "finished constructing weighting factors"
            << endl;
    }
}


void solidPolyMesh::calcCellImagEdges() const
{
    const polyMesh& mesh = (*this)();

    cellImagEdgesPtr_ = new labelListList(mesh.nCells());
    labelListList& cellImagEdges = *cellImagEdgesPtr_;

    const labelListList& cpc = cellPointCells();
    labelList imEdgesFound(mesh.nCells(), 0);

    forAll(cellImagEdges, cellI)
    {
        cellImagEdges[cellI].setSize(cpc[cellI].size(), -1);
    }

    const edgeList& imEdges = imagEdges();

    forAll(imEdges, imEdgeI)
    {
        const edge& e = imEdges[imEdgeI];
        const label startCellID = e.start();
        const label endCellID = e.end();

        cellImagEdges[startCellID][imEdgesFound[startCellID]] = imEdgeI;
        cellImagEdges[endCellID][imEdgesFound[endCellID]] = imEdgeI;

        imEdgesFound[startCellID]++;
        imEdgesFound[endCellID]++;
    }

    // Check all edges were set
    label minImEdgeID = 0;
    forAll(cellImagEdges, cellI)
    {
        minImEdgeID = min(cellImagEdges[cellI]);
    }

    if (minImEdgeID < 0)
    {
        FatalErrorIn
        (
            "void solidPolyMesh::calcCellImagEdges() const"
        )   << "Not all cell imaginary edges set" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
solidPolyMesh::solidPolyMesh(const polyMesh& pMesh)
:
    GeoMesh<polyMesh>(pMesh),
    boundary_(*this, pMesh.boundaryMesh()),
    lduPtr_(NULL),
    imagEdgesPtr_(NULL),
    nearImagEdgesPtr_(NULL),
    fvMeshAddressingMapPtr_(NULL),
    pointWeightsPtr_(NULL),
    edgeWeightsPtr_(NULL),
    cellPointCellsPtr_(NULL),
    cellImagEdgesPtr_(NULL)
{
    if (debug)
    {
        Info<< "solidPolyMesh::solidPolyMesh(const polyMesh&) : "
            << "Creating solidPolyMesh" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

solidPolyMesh::~solidPolyMesh()
{
    if (debug)
    {
        Info<< "solidPolyMesh::~solidPolyMesh() : "
            << "Deleting solidPolyMesh" << endl;
    }

    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const lduAddressing& solidPolyMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new solidPolyMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


const edgeList& solidPolyMesh::imagEdges() const
{
    if (!imagEdgesPtr_)
    {
        calcImagEdges();
    }

    return *imagEdgesPtr_;
}


const boolList& solidPolyMesh::nearImagEdges() const
{
    if (!nearImagEdgesPtr_)
    {
        calcImagEdges();
    }

    return *nearImagEdgesPtr_;
}


const labelList& solidPolyMesh::fvMeshAddressingMap() const
{
    if (!fvMeshAddressingMapPtr_)
    {
        calcFvMeshAddressingMap();
    }

    return *fvMeshAddressingMapPtr_;
}


const scalarListList& solidPolyMesh::pointWeights() const
{
    if (!pointWeightsPtr_)
    {
        makePointWeights();
    }

    return *pointWeightsPtr_;
}


const scalarListList& solidPolyMesh::edgeWeights() const
{
    if (!edgeWeightsPtr_)
    {
        makeEdgeWeights();
    }

    return *edgeWeightsPtr_;
}


const labelListList& solidPolyMesh::cellPointCells() const
{
    if (!cellPointCellsPtr_)
    {
        calcImagEdges();
    }

    return *cellPointCellsPtr_;
}


const labelListList& solidPolyMesh::cellImagEdges() const
{
    if (!cellImagEdgesPtr_)
    {
        calcCellImagEdges();
    }

    return *cellImagEdgesPtr_;
}


label solidPolyMesh::findImaginaryEdge
(
    const label firstCellID,
    const label secondCellID
) const
{
    // Check which cell has fewer imag edges connected to it
    label curCellID = firstCellID;
    label otherCellID = secondCellID;
    if
    (
        cellImagEdges()[firstCellID].size()
        > cellImagEdges()[secondCellID].size()
    )
    {
        curCellID = secondCellID;
        otherCellID = firstCellID;
    }
    const labelList& curCellImagEdges = cellImagEdges()[curCellID];

    // We will loop through all imaginary edges of curCell until we find the
    // edge which connects the curCell and otherCell

    label foundEdgeID = -1;
    const edgeList& imEdges = imagEdges();

    forAll(curCellImagEdges, imEdgeI)
    {
        const edge curImEdge = imEdges[curCellImagEdges[imEdgeI]];
        if
        (
            curImEdge.start() == otherCellID || curImEdge.end() == otherCellID
        )
        {
            foundEdgeID = imEdgeI;
            break;
        }
    }

    if (foundEdgeID == -1)
    {
        return -1;
    }

    return curCellImagEdges[foundEdgeID];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool solidPolyMesh::operator!=(const solidPolyMesh& bm) const
{
    return &bm != this;
}


bool solidPolyMesh::operator==(const solidPolyMesh& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
