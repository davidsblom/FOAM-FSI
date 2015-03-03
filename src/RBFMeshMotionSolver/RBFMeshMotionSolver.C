
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFMeshMotionSolver.H"

using namespace Foam;

defineTypeNameAndDebug( RBFMeshMotionSolver, 0 );

addToRunTimeSelectionTable
(
  motionSolver,
  RBFMeshMotionSolver,
  dictionary
);

RBFMeshMotionSolver::RBFMeshMotionSolver(
  const polyMesh & mesh,
  Istream & msData
  )
  :
  motionSolver( mesh ),
  motionCenters( mesh.boundaryMesh().size(), vectorField( 0 ) ),
  staticPatches( lookup( "staticPatches" ) ),
  staticPatchIDs( staticPatches.size() ),
  movingPatches( lookup( "movingPatches" ) ),
  movingPatchIDs( movingPatches.size() ),
  newPoints( mesh.points().size(), vector::zero ),
  nbGlobalFaceCenters( Pstream::nProcs(), 0 ),
  nbGlobalMovingFaceCenters( Pstream::nProcs(), 0 ),
  nbGlobalStaticFaceCenters( Pstream::nProcs(), 0 ),
  twoDCorrector( mesh ),
  nbPoints( 0 )
{
  // Find IDs of staticPatches_
  forAll( staticPatches, patchI )
  {
    label patchIndex = mesh.boundaryMesh().findPatchID( staticPatches[patchI] );

    assert( patchIndex >= 0 );

    staticPatchIDs[patchI] = patchIndex;
  }

  // Find IDs of movingPatches_
  forAll( movingPatches, patchI )
  {
    label patchIndex = mesh.boundaryMesh().findPatchID( movingPatches[patchI] );

    assert( patchIndex >= 0 );

    movingPatchIDs[patchI] = patchIndex;
  }

  // Verify that a patch is not defined as a static and a moving patch

  forAll( staticPatchIDs, staticPatchI )
  {
    // Search the moving patches for static patchI
    forAll( movingPatchIDs, movingPatchI )
    {
      assert( movingPatchIDs[movingPatchI] != staticPatchIDs[staticPatchI] );
    }
  }

  // Initialize RBF interpolator

  dictionary & dict = subDict( "interpolation" );

  word function = dict.lookup( "function" );

  assert( function == "TPS" || function == "WendlandC0" || function == "WendlandC2" || function == "WendlandC4" || function == "WendlandC6" );

  std::shared_ptr<rbf::RBFFunctionInterface> rbfFunction;

  Info << "Radial Basis Function interpolation: Selecting RBF function: " << function << endl;

  if ( function == "TPS" )
    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::TPSFunction() );

  if ( function == "WendlandC0" )
  {
    scalar radius = readScalar( dict.lookup( "radius" ) );
    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC0Function( radius ) );
  }

  if ( function == "WendlandC2" )
  {
    scalar radius = readScalar( dict.lookup( "radius" ) );
    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC2Function( radius ) );
  }

  if ( function == "WendlandC4" )
  {
    scalar radius = readScalar( dict.lookup( "radius" ) );
    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC4Function( radius ) );
  }

  if ( function == "WendlandC6" )
  {
    scalar radius = readScalar( dict.lookup( "radius" ) );
    rbfFunction = std::shared_ptr<rbf::RBFFunctionInterface> ( new rbf::WendlandC6Function( radius ) );
  }

  assert( rbfFunction );

  std::shared_ptr<rbf::RBFInterpolation> rbfInterpolator( new rbf::RBFInterpolation( rbfFunction ) );

  bool coarsening = readBool( subDict( "coarsening" ).lookup( "enabled" ) );
  double tol = readScalar( subDict( "coarsening" ).lookup( "tol" ) );
  bool livePointSelection = readBool( subDict( "coarsening" ).lookup( "livePointSelection" ) );
  double tolLivePointSelection = readScalar( subDict( "coarsening" ).lookup( "tolLivePointSelection" ) );
  int coarseningMinPoints = readLabel( subDict( "coarsening" ).lookup( "minPoints" ) );
  int coarseningMaxPoints = readLabel( subDict( "coarsening" ).lookup( "maxPoints" ) );
  bool exportSelectedPoints = readBool( subDict( "coarsening" ).lookup( "exportSelectedPoints" ) );
  rbf = std::shared_ptr<rbf::RBFCoarsening> ( new rbf::RBFCoarsening( rbfInterpolator, coarsening, livePointSelection, true, tol, tolLivePointSelection, coarseningMinPoints, coarseningMaxPoints, exportSelectedPoints ) );
}

RBFMeshMotionSolver::~RBFMeshMotionSolver()
{}

tmp<pointField> RBFMeshMotionSolver::curPoints() const
{
  // Prepare new points: same as old point
  tmp<pointField> tnewPoints
  (
    new vectorField( mesh().nPoints(), vector::zero )
  );

  pointField & newPoints = tnewPoints();

  newPoints = this->newPoints;

  // Add old point positions
  newPoints += mesh().points();

  return tnewPoints;
}

// As a first step, the motion is defined in the
void RBFMeshMotionSolver::setMotion( const Field<vectorField> & motion )
{
  // Input checking

  assert( motion.size() == mesh().boundaryMesh().size() );

  forAll( motion, ipatch )
  {
    const vectorField & mpatch = motion[ipatch];

    // Check whether the size of patch motion is equal to number of face centers in patch
    if ( mpatch.size() > 0 )
      assert( mpatch.size() == mesh().boundaryMesh()[ipatch].faceCentres().size() );

    // Check whether the size of a moving patch is equal to the number of face centers in the patch
    // First check if patchid is a moving patch
    bool movingPatch = false;
    forAll( movingPatchIDs, movingPatchI )
    {
      if ( movingPatchIDs[movingPatchI] == ipatch )
        movingPatch = true;
    }

    if ( movingPatch )
      assert( mpatch.size() == mesh().boundaryMesh()[ipatch].faceCentres().size() );
  }

  motionCenters = motion;
}

void RBFMeshMotionSolver::updateMesh( const mapPolyMesh & )
{
  assert( false );
}

void RBFMeshMotionSolver::solve()
{
  assert( motionCenters.size() == mesh().boundaryMesh().size() );

  /*
   * RBF interpolator from face centers to local complete mesh vertices
   * The interpolation consists of the following steps:
   * 1. Build a matrix with the face center positions of the static patches and the moving patches
   * 2. Build a matrix with the positions of every vertex in the local mesh
   * 3. Build a matrix with the displacement/motion of the face center positions of the static patches and the moving patches
   * 4. Perform the interpolation from the face centers to the complete mesh
   * 5. Correct the mesh vertices of the static patches. Set these displacement to zero.
   * 6. Set the motion of the mesh vertices
   */

  /*
   * Step 1: Build a matrix with the face center positions of the static patches and the moving patches
   * The order of the matrix is defined as first a list of the moving patch face centers,
   * thereafter the static patch face centers. These are the control points used by the
   * radial basis function interpolation.
   * The control points should be exactly the same at each processor, and are therefore communicated
   * to each process. As only an absolute RBF interpolation is implemented, this communication step is only
   * performed once per simulation.
   * The global ordering of the data is first the information of the moving patches,
   * thereafter the static patches.
   */

  int nbFaceCenters = 0;
  int nbMovingFaceCenters = 0;
  int nbStaticFaceCenters = 0;

  if ( sum( nbGlobalFaceCenters ) == 0 )
  {
    // Determine the number of face centers
    // The total number of face centers is simply the sum of the face centers
    // on each processor.
    nbFaceCenters = 0;

    // First add the moving patches, and thereafter the static patches

    forAll( movingPatchIDs, i )
    {
      nbMovingFaceCenters += mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres().size();
    }

    forAll( staticPatchIDs, i )
    {
      nbStaticFaceCenters += mesh().boundaryMesh()[staticPatchIDs[i]].faceCentres().size();
    }

    // Calculate sum of all faces on each processor
    nbGlobalMovingFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters;
    nbGlobalStaticFaceCenters[Pstream::myProcNo()] = nbStaticFaceCenters;
    nbGlobalFaceCenters[Pstream::myProcNo()] = nbMovingFaceCenters + nbStaticFaceCenters;

    reduce( nbGlobalMovingFaceCenters, sumOp<labelList>() );
    reduce( nbGlobalStaticFaceCenters, sumOp<labelList>() );
    reduce( nbGlobalFaceCenters, sumOp<labelList>() );
  }

  nbMovingFaceCenters = sum( nbGlobalMovingFaceCenters );
  nbStaticFaceCenters = sum( nbGlobalStaticFaceCenters );
  nbFaceCenters = sum( nbGlobalFaceCenters );

  // Determine the offset taking into account multiple processors

  int globalMovingOffset = 0;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    globalMovingOffset += nbGlobalMovingFaceCenters[i];

  int globalStaticOffset = nbMovingFaceCenters;

  for ( int i = 0; i < Pstream::myProcNo(); i++ )
    globalStaticOffset += nbGlobalStaticFaceCenters[i];

  if ( !rbf->rbf->computed )
  {
    rbf::matrix positions( nbFaceCenters, mesh().nGeometricD() );
    positions.setZero();
    int offset = 0;

    vectorField positionsField( positions.rows(), vector::zero );

    forAll( movingPatchIDs, i )
    {
      const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

      // Set the positions for patch i
      forAll( faceCentres, j )
      {
        positionsField[j + offset + globalMovingOffset] = faceCentres[j];
      }

      offset += faceCentres.size();
    }

    offset = 0;

    forAll( staticPatchIDs, i )
    {
      const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[staticPatchIDs[i]].faceCentres();

      // Set the positions for patch i
      forAll( faceCentres, j )
      {
        positionsField[j + offset + globalStaticOffset] = faceCentres[j];
      }

      offset += faceCentres.size();
    }

    reduce( positionsField, sumOp<vectorField>() );

    // Copy the FOAM vector field to an Eigen matrix
    for ( int i = 0; i < positions.rows(); i++ )
      for ( int j = 0; j < positions.cols(); j++ )
        positions( i, j ) = positionsField[i][j];

    /*
     * Step 2: Build a matrix with the positions of every vertex in the local mesh.
     * This is only local information and does not need to be communicated to other
     * processors.
     */

    const Foam::pointField & points = mesh().points();

    // Determine the number of points by using the 2d corrector
    nbPoints = 0;
    forAll( points, i )
    {
      if ( twoDCorrector.marker()[i] == 0 )
        nbPoints++;
    }

    rbf::matrix positionsInterpolation( nbPoints, positions.cols() );

    int index = 0;
    forAll( points, i )
    {
      if ( twoDCorrector.marker()[i] == 0 )
      {
        for ( int j = 0; j < positionsInterpolation.cols(); j++ )
          positionsInterpolation( index, j ) = points[i][j];

        index++;
      }
    }

    // rbf->compute( positions, positionsInterpolation, positionsParallelLocation );
    rbf->compute( positions, positionsInterpolation );

    rbf->setNbMovingAndStaticFaceCenters( nbMovingFaceCenters, nbStaticFaceCenters );
  }

  /*
   * Step 3: Build a matrix with the displacement/motion of the face center
   * positions of the static patches and the moving patches.
   * The motion needs to be communicated to every process at every mesh deformation.
   * This is considered to be the most expensive step with regards to parallel
   * scalability of the overall algorithm.
   */

  rbf::matrix values( nbFaceCenters, mesh().nGeometricD() );
  values.setZero();

  vectorField valuesField( values.rows(), vector::zero );

  int offset = 0;

  forAll( movingPatchIDs, i )
  {
    const Foam::vectorField::subField faceCentres = mesh().boundaryMesh()[movingPatchIDs[i]].faceCentres();

    forAll( motionCenters[movingPatchIDs[i]], j )
    {
      valuesField[j + offset + globalMovingOffset] = motionCenters[movingPatchIDs[i]][j];
    }

    offset += faceCentres.size();
  }

  reduce( valuesField, sumOp<vectorField>() );

  // Copy the FOAM vector field to an Eigen matrix
  for ( int i = 0; i < values.rows(); i++ )
    for ( int j = 0; j < values.cols(); j++ )
      values( i, j ) = valuesField[i][j];

  /*
   * Step 4: Perform the interpolation from the face centers to the complete mesh
   */

  rbf::matrix valuesInterpolation( nbPoints, values.cols() );
  valuesInterpolation.setZero();

  rbf->interpolate( values, valuesInterpolation );

  // Apply the 2d correction

  vectorField valuesInterpolationField( mesh().points().size(), Foam::vector::zero );
  int index = 0;
  forAll( valuesInterpolationField, i )
  {
    if ( twoDCorrector.marker()[i] == 0 )
    {
      for ( int j = 0; j < valuesInterpolation.cols(); j++ )
        valuesInterpolationField[i][j] = valuesInterpolation( index, j );

      index++;
    }
  }

  twoDCorrector.setShadowSide( valuesInterpolationField );

  /*
   * Step 5: Correct the mesh vertices of the static patches. Set these displacements to zero.
   */

  // Loop over all the patches, and set the static patches to zero.
  forAll( mesh().boundaryMesh(), i )
  {
    const labelList & meshPoints = mesh().boundaryMesh()[i].meshPoints();

    bool isStaticPatch = false;
    forAll( staticPatchIDs, j )
    {
      if ( i == staticPatchIDs[j] )
        isStaticPatch = true;
    }

    if ( isStaticPatch )
    {
      for ( int j = 0; j < meshPoints.size(); j++ )
        valuesInterpolationField[meshPoints[j]] = Foam::vector::zero;
    }

    offset += meshPoints.size();
  }

  /*
   * Step 6: Set the motion of the mesh vertices
   */

  assert( newPoints.size() == valuesInterpolationField.size() );

  newPoints = valuesInterpolationField;
}
