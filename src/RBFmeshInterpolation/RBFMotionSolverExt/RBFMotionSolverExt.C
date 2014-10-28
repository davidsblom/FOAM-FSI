
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFMotionSolverExt.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "SortableList.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug( RBFMotionSolverExt, 0 );

  addToRunTimeSelectionTable
  (
    motionSolver,
    RBFMotionSolverExt,
    dictionary
  );

  // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

  void RBFMotionSolverExt::makeControlIDs()
  {
    const pointField & points = mesh().points();
    labelList markedPoints( points.size(), 0 );
    internalIDs_.setSize( points.size() );
    internalPoints_.setSize( points.size() );

    // If parallel run make sure that only the owner of processor patches get the static boundary points
    if ( Pstream::parRun() )
    {
      const polyBoundaryMesh & bMesh( mesh().boundaryMesh() );
      const labelList & procPatches( mesh().globalData().processorPatches() );
      forAll( procPatches, i )       // For all processor patches
      {
        const processorPolyPatch & procPatch = refCast<const processorPolyPatch>( bMesh[procPatches[i]] );
        const labelList & procPatchPoints = procPatch.meshPoints();

        if ( procPatch.neighbour() )             // If not the owner of the patch
        {
          forAll( procPatchPoints, k )                   // Go over all points and mark them with -5
          {
            markedPoints[procPatchPoints[k]] = -5;
          }
        }
      }
    }

    // =======================================================================//
    // ============================= Mark points =============================//
    // =======================================================================//

    label nMarkedPoints = 0;

    // Start with movingpatches to ensure that if point is shared by moving and static patch it is marked as static
    // * * * * * * * * * mark moving points * * * * * * * * * //
    forAll( movingPatches_, patchI )
    {
      label patchIndex = movingPatchIDs_[patchI];
      const labelList & mp = mesh().boundaryMesh()[patchIndex].meshPoints();

      forAll( mp, i )
      {
        // Check if on neigbour side of proc patch (-5 or already 2)
        if ( markedPoints[mp[i]] == -5 || markedPoints[mp[i]] == 2 )
        {
          markedPoints[mp[i]] = 2;              // Moving points but also on neighbour side of proc patch
        }
        else
        {
          markedPoints[mp[i]] = 1;              // Moving point
        }

        nMarkedPoints++;
      }
    }

    // Initialize sizes of relevant variables
    movingIDs_.setSize( nMarkedPoints );
    movingPoints_.setSize( nMarkedPoints );

    // * * * * * * * * * mark static points * * * * * * * * * //
    nMarkedPoints = 0;     // reset counter
    forAll( staticPatches_, patchI )
    {
      label patchIndex = staticPatchIDs_[patchI];
      const labelList & mp = mesh().boundaryMesh()[patchIndex].meshPoints();

      forAll( mp, i )
      {
        // Check if on neighbour side of proc patch (both moving marked (2) and just neighbour proc patch marked (-5 || -2))
        if ( markedPoints[mp[i]] == -5 || markedPoints[mp[i]] == 2 || markedPoints[mp[i]] == -2 )
        {
          markedPoints[mp[i]] = -2;
        }
        else
        {
          markedPoints[mp[i]] = -1;              // This also overrides any moving points identified before (static patch points overrule moving patch points)
        }

        nMarkedPoints++;
      }
    }

    // Initialize sizes of relevant variables
    staticIDs_.setSize( nMarkedPoints );
    staticPoints_.setSize( nMarkedPoints );
    controlIDs_.setSize( movingIDs_.size() + staticIDs_.size() );
    controlPoints_.setSize( movingIDs_.size() + staticIDs_.size() );

    // =======================================================================//
    // ============================= Point IDs and points ====================//
    // =======================================================================//
    label nMovingPoints = 0;
    label nStaticPoints = 0;
    label nInternalPoints = 0;

    // * * * * * * * * * set moving IDs & points * * * * * * * * * //
    // In movingIDs_ and movingPoints_ duplicates might be available. However, at picking the control points this is taken into account
    // For moving points this is also more relax, since setMotion(m) is based on the patch.localPoints().
    // This can cause duplicates when mutiple connected moving patches are used. In this way the motion is enfored properly without
    // destroying the correct control point selectionpatchMovingGamma_
    forAll( movingPatchIDs_, patchI )
    {
      label patchIndex = movingPatchIDs_[patchI];
      const labelList & mp = mesh().boundaryMesh()[patchIndex].meshPoints();
      patchMovingGamma_[patchIndex] = labelList( mp.size(), 0 );
      forAll( mp, i )
      {
        // If moving: marker=1 or =2 (on neighbour side of proc patch, but still moving point->never control point)
        if ( markedPoints[mp[i]] == 1 || markedPoints[mp[i]] == 2 )
        {
          // Set moving ID and point
          movingIDs_[nMovingPoints] = mp[i];
          movingPoints_[nMovingPoints] = points[mp[i]];      // Not in original->do not remember why
          patchMovingGamma_[patchIndex][i] = 1;
          nMovingPoints++;
        }
      }
    }

    // * * * * * * * * * set static IDs & points * * * * * * * * * //
    // In staticIDs_ and staticPoints_ duplicates might be available. However, at picking the control points this is taken into account
    forAll( staticPatchIDs_, patchI )
    {
      label patchIndex = staticPatchIDs_[patchI];
      const labelList & mp = mesh().boundaryMesh()[patchIndex].meshPoints();

      forAll( mp, i )
      {
        // If static: marker=-1 or =-2 (on neighbour side of proc patch, but still static point->never control point)
        if ( markedPoints[mp[i]] == -1 || markedPoints[mp[i]] == -2 )
        {
          // Set static ID and points
          staticIDs_[nStaticPoints] = mp[i];
          staticPoints_[nStaticPoints] = points[mp[i]];
          nStaticPoints++;
        }
      }
    }

    // * * * * * * * * * set internal IDs & points * * * * * * * * * //
    forAll( markedPoints, i )
    {
      // If internal: marker=0 or -5 (not marked or marked as neighbour side proc patch, but not on either static or moving patches)
      if ( (markedPoints[i] == 0 || markedPoints[i] == -5) && twoDPointCorrectorRBF_.marker()[i] == 0 )
      {
        // if(markedPoints[i] == 0 || markedPoints[i] == -5)

        // Set internal ID and points
        internalIDs_[nInternalPoints] = i;
        internalPoints_[nInternalPoints] = points[i];
        markedPoints[i] = 0;        // Reset to zero. If -5 it is not on moving or static patch so can be set back to 0
        nInternalPoints++;
      }
    }

    // Set sizes accordingly
    movingIDs_.setSize( nMovingPoints );
    movingPoints_.setSize( nMovingPoints );
    staticIDs_.setSize( nStaticPoints );
    staticPoints_.setSize( nStaticPoints );
    internalIDs_.setSize( nInternalPoints );
    internalPoints_.setSize( nInternalPoints );

    // =======================================================================//
    // ============================= Control IDs and points ====================//
    // =======================================================================//

    // * * * * * * * * * pick static control points * * * * * * * * * //
    label nStaticControlPoints = 0; // Counts all static control points
    forAll( staticPatches_, patchI )
    {
      label nPickedPoints = 0;  // Count picked points of this patch
      label patchIndex = staticPatchIDs_[patchI];
      const labelList & mp = mesh().boundaryMesh()[patchIndex].meshPoints();

      label possiblePoints = 0;
      forAll( mp, i )
      {
        // This is for twoD cases-> only pick from one plane
        if ( twoDPointCorrectorRBF_.marker()[mp[i]] == 0 )
        {
          // Pick only when marker == -1, which means it is a static point
          if ( markedPoints[mp[i]] == -1 )
          {
            // Apply Coarsening
            if ( possiblePoints >= nPickedPoints * coarseningRatioStat_ )
            {
              controlIDs_[nStaticControlPoints] = mp[i];
              controlPoints_[nStaticControlPoints] = points[mp[i]];
              markedPoints[mp[i]] = -3;                  // prevents points to be picked twice if they live on multiple (static) patches

              // Update counters
              nPickedPoints++;
              nStaticControlPoints++;
            }

            possiblePoints++;
          }
        }
      }
    }
    nStaticControlPoints_ = nStaticControlPoints;

    // * * * * * * * * * pick moving control points * * * * * * * * * //
    label nMovingControlPoints = 0;    // Counts all moving control points
    forAll( movingPatches_, patchI )
    {
      label nPickedPoints = 0;  // Count picked points of this patch
      label patchIndex = movingPatchIDs_[patchI];
      const labelList & mp = mesh().boundaryMesh()[patchIndex].meshPoints();

      label possiblePoints = 0;
      forAll( mp, i )
      {
        // This is for twoD cases-> only pick from one plane
        if ( twoDPointCorrectorRBF_.marker()[mp[i]] == 0 )
        {
          // Pick only when marker == 3, which means it is a moving point
          if ( markedPoints[mp[i]] == 1 )
          {
            // /Apply Coarsening
            if ( possiblePoints >= nPickedPoints * coarseningRatioMov_ )
            {
              controlIDs_[nStaticControlPoints + nMovingControlPoints] = mp[i];
              controlPoints_[nStaticControlPoints + nMovingControlPoints] = points[mp[i]];
              markedPoints[mp[i]] = 3;          // prevents points to be picked twice if they live on multiple (moving) patches

              // Update counters
              nPickedPoints++;
              nMovingControlPoints++;
            }

            possiblePoints++;
          }
        }
      }
    }

    // Set sizes and summarize
    controlIDs_.setSize( nStaticControlPoints + nMovingControlPoints );
    controlPoints_.setSize( nStaticControlPoints + nMovingControlPoints );

    // ======================================================================================//
    // ============================= Internal Points ========================================//
    // ======================================================================================//

    nInternalPoints = 0; // reuse counter
    forAll( markedPoints, i )
    {
      if ( markedPoints[i] == 0 && twoDPointCorrectorRBF_.marker()[i] == 0 ) // Only if zero since -4 should not exist anymore (see internalIDs filling)
      {
        internalIDs_[nInternalPoints] = i;
        internalPoints_[nInternalPoints] = points[i];
        nInternalPoints++;
      }
    }
    internalIDs_.setSize( nInternalPoints );
    internalPoints_.setSize( nInternalPoints );
  }

  void RBFMotionSolverExt::moveMesh( pointField & newPoints ) const
  {
    vectorField movingControlPointsMotion( sum( nGlobalMovingControlPoints_ ), vector::zero );

    // Calculate starting index in movingControlPointsMotion for each proc
    labelList startInd( Pstream::nProcs(), 0 );

    for ( int i = 1; i < Pstream::nProcs(); i++ )
    {
      startInd[i] = startInd[i - 1] + nGlobalMovingControlPoints_[i - 1];
    }

    // Fill vectorField movingControlPointsMotion on each proc at the right location (also works for serial calculations)
    for ( int i = 0; i < nGlobalMovingControlPoints_[Pstream::myProcNo()]; i++ )
    {
      movingControlPointsMotion[startInd[Pstream::myProcNo()] + i] = newPoints[movingControlIDs_[i]];
    }

    // And reduce it to a single vectorField on all procs
    reduce( movingControlPointsMotion, sumOp<vectorField>() );

    // Call interpolation
    vectorField interpolatedMotion = interpolationRed_.interpolateRed( movingControlPointsMotion );

    // Set interpolated values
    forAll( internalIDs_, i )
    {
      newPoints[internalIDs_[i]] = interpolatedMotion[i];
    }
  }

  void RBFMotionSolverExt::setAndOrderPatchIDs()
  {
    // Find IDs of staticPatches_
    forAll( staticPatches_, patchI )
    {
      label patchIndex = mesh().boundaryMesh().findPatchID( staticPatches_[patchI] );

      if ( patchIndex < 0 )
      {
        FatalErrorIn( "void RBFMotionSolver::setAndOrderPatchIDs()" )
        << "Patch " << staticPatches_[patchI] << " not found.  "
        << "valid patch names: " << mesh().boundaryMesh().names()
        << abort( FatalError );
      }
      else
      {
        staticPatchIDs_[patchI] = patchIndex;
      }
    }

    // Find IDs of movingPatches_
    forAll( movingPatches_, patchI )
    {
      label patchIndex = mesh().boundaryMesh().findPatchID( movingPatches_[patchI] );

      if ( patchIndex < 0 )
      {
        FatalErrorIn( "void RBFMotionSolver::setAndOrderPatchIDs()" )
        << "Patch " << movingPatches_[patchI] << " not found.  "
        << "valid patch names: " << mesh().boundaryMesh().names()
        << abort( FatalError );
      }
      else
      {
        movingPatchIDs_[patchI] = patchIndex;
      }
    }

    // Sort movingPatches_ on patchIDs
    SortableList<label> sortedMovingPatchIDs( movingPatchIDs_.size() );
    sortedMovingPatchIDs = movingPatchIDs_;
    sortedMovingPatchIDs.sort();

    // Order movingPatches accordingly
    wordList tMovingPatches = movingPatches_;
    forAll( movingPatchIDs_, i )
    {
      tMovingPatches[i] = movingPatches_[sortedMovingPatchIDs.indices()[i]];
    }

    // Copy data back to holders
    movingPatches_ = tMovingPatches;
    movingPatchIDs_ = sortedMovingPatchIDs;
  }

  // Get back the point IDs from the moving control points for each PROC
  void RBFMotionSolverExt::setMovingControlPointsIDs()
  {
    label nMovingControlPoints = controlIDs_.size() - nStaticControlPoints_;

    movingControlIDs_.setSize( nMovingControlPoints );

    // Set control IDs in right order (patchID order and patch.meshpoint order) for each PROC
    for ( label i = nStaticControlPoints_; i < controlIDs_.size(); i++ )
    {
      // Moving control points are the last control points
      movingControlIDs_[i - nStaticControlPoints_] = controlIDs_[i];
    }
  }

  // Get back the indices of the moving control points in the controlPoints_ list (for both serial and par)
  const labelList RBFMotionSolverExt::movingControlPointsIndices() const
  {
    const label nGlobalStaticControlPoints = sum( nGlobalStaticControlPoints_ );  // Global number of static control points
    const label nGlobalMovingControlPoints = sum( nGlobalMovingControlPoints_ );

    labelList movingControlPointsIndex( nGlobalMovingControlPoints );

    for ( label i = nGlobalStaticControlPoints; i < controlPoints_.size(); i++ )
    {
      // Last controlPoints are the moving control points
      movingControlPointsIndex[i - nGlobalStaticControlPoints] = i;
    }

    return movingControlPointsIndex;
  }

  void RBFMotionSolverExt::gatherControlPointsExt()
  {
    // Get number of control points on each proc (total and moving)
    labelList nControlPoints( Pstream::nProcs(), 0 );
    nControlPoints[Pstream::myProcNo()] = controlIDs_.size();

    // Reduce (gather, scatter) over procs
    reduce( nControlPoints, sumOp<labelList>() );
    reduce( nGlobalStaticControlPoints_, sumOp<labelList>() );
    reduce( nGlobalMovingControlPoints_, sumOp<labelList>() );

    // Get some values for further creation of global controlPoint list
    label nTotalControlPoints = sum( nControlPoints );
    label nTotalMovingControlPoints = sum( nGlobalMovingControlPoints_ );

    // Determine starting index for each proc
    labelList startingIndStat( Pstream::nProcs(), 0 );
    labelList startingIndMov( Pstream::nProcs(), nTotalControlPoints - nTotalMovingControlPoints );

    for ( int i = 1; i < Pstream::nProcs(); i++ )
    {
      startingIndStat[i] = startingIndStat[i - 1] + (nControlPoints[i - 1] - nGlobalMovingControlPoints_[i - 1]);
      startingIndMov[i] = startingIndMov[i - 1] + nGlobalMovingControlPoints_[i - 1];
    }

    // Set control points at right lcoation in globalControlPoint vector
    vectorField globalControlPoints( nTotalControlPoints, vector::zero );

    // Static control point in first part
    for ( int i = 0; i < nStaticControlPoints_; i++ )
    {
      label ind = startingIndStat[Pstream::myProcNo()] + i;
      globalControlPoints[ind] = controlPoints_[i];
    }

    // Moving control points in second part
    for ( int i = 0; i < controlIDs_.size() - nStaticControlPoints_; i++ )
    {
      label ind = startingIndMov[Pstream::myProcNo()] + i;
      globalControlPoints[ind] = controlPoints_[nStaticControlPoints_ + i];
    }

    // Ensure all procs has the same globalControlPoint list
    reduce( globalControlPoints, sumOp<vectorField>() );
    controlPoints_ = globalControlPoints;    // Set controlPoints_ at all procs to globalControlPoints
  }

  void RBFMotionSolverExt::prepareParallelRun()
  {
    if ( Pstream::parRun() )
    {
      gatherControlPointsExt();
    }
  }

  void RBFMotionSolverExt::checkInterpolationDimension()
  {
    word dim( subDict( "interpolation" ).lookup( "Dim" ) );

    if ( twoDPointCorrectorRBF_.required() )
    {
      if ( dim == "ThreeD" )
      {
        FatalErrorIn( "void RBFMotionSolverExt::setInterpolationDimension" )
        << "Incorrect dimension according to mesh: mesh = TwoD, while interpolation = ThreeD" << endl
        << "Please correct this." << endl
        << abort( FatalError );
      }
    }
    else
    {
      if ( dim == "TwoD" )
      {
        FatalErrorIn( "void RBFMotionSolverExt::setInterpolationDimension" )
        << "Incorrect dimension according to mesh: mesh = ThreeD, while interpolation = TwoD" << endl
        << "Please correct this." << endl
        << abort( FatalError );
      }
    }
  }

  // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

  RBFMotionSolverExt::RBFMotionSolverExt(
    const polyMesh & mesh,
    Istream &
    )
    :
    motionSolver( mesh ),
    staticPatches_( lookup( "staticPatches" ) ),
    staticPatchIDs_( staticPatches_.size() ),
    movingPatches_( lookup( "movingPatches" ) ),
    movingPatchIDs_( movingPatches_.size() ),
    coarseningRatioStat_( readLabel( lookup( "coarseningRatioStat" ) ) ),
    coarseningRatioMov_( readLabel( lookup( "coarseningRatioMov" ) ) ),
    controlIDs_( 0 ),
    controlPoints_( 0 ),
    internalIDs_( 0 ),
    internalPoints_( 0 ),
    staticIDs_( 0 ),
    staticPoints_( 0 ),
    movingIDs_( 0 ),
    movingPoints_( 0 ),
    patchMovingGamma_( mesh.boundaryMesh().size(), labelField( 0 ) ),
    motion_( 0 ),
    movingControlPointsIndex_( 0 ),
    movingControlIDs_( 0 ),
    movingPointsPatchStart_( movingPatches_.size() ),
    nStaticControlPoints_( 0 ),
    nGlobalMovingControlPoints_( Pstream::nProcs(), 0 ),
    nGlobalStaticControlPoints_( Pstream::nProcs(), 0 ),
    interpolationRed_
    (
    mesh,
    subDict( "interpolation" ),
    controlPoints_,
    internalPoints_
    ),
    twoDPointCorrectorRBF_( mesh )
  {
    // Set patch IDs of moving and static patches and order moving patches in ID order (low to high)
    setAndOrderPatchIDs();

    // Check whether interpolation dimensions match with mesh dimensions
    checkInterpolationDimension();

    makeControlIDs();    // Set IDs for all (local) points
    setMovingControlPointsIDs();    // Set movingControlIDs_ based on the controlIDs and nStaticControlPoints_

    // Used to gather all control points properly
    nGlobalStaticControlPoints_[Pstream::myProcNo()] = nStaticControlPoints_;
    nGlobalMovingControlPoints_[Pstream::myProcNo()] = controlIDs_.size() - nStaticControlPoints_;

    // Gathers the controlpoints and scatters them on all procs
    prepareParallelRun();

    motion_.setSize( movingIDs_.size(), vector::zero );

    // Construct matrix for RBF interpolation by indicating which controlPoints_ are moving
    labelList movingControlPointsIndex = movingControlPointsIndices();
    interpolationRed_.createReducedEvaluationMatrix( movingControlPointsIndex );
  }

  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  RBFMotionSolverExt::~RBFMotionSolverExt()
  {}


  // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void RBFMotionSolverExt::setMotion( const Field<vectorField> & m )
  {
    if ( m.size() != mesh().boundaryMesh().size() )
    {
      FatalErrorIn( "void RBFMotionSolverExt::setMotion(const Field<vectorField>& m)" )
      << "Incorrect size of motion Field: m = " << m.size()
      << " Number of patches = " << mesh().boundaryMesh().size()
      << abort( FatalError );
    }

    label motionCounter = 0;
    forAll( m, ipatch )
    {
      const vectorField & mpatch = m[ipatch];

      // Check whether the size of patch motion is equal to number of points in patch
      if ( mpatch.size() > 0 && mpatch.size() != mesh().boundaryMesh()[ipatch].meshPoints().size() )
      {
        FatalErrorIn( "void RBFMotionSolverExt::setMotion(const Field<vectorField>& m)" )
        << "Incorrect size of motion Field: m[" << ipatch << "] = " << m[ipatch].size()
        << " Number of point ins patch = " << mesh().boundaryMesh()[ipatch].meshPoints().size()
        << abort( FatalError );
      }

      const labelField & gamma = patchMovingGamma_[ipatch];

      // Loop over patch points and only assign the motion when it is identified as motionPoint in this class
      forAll( mpatch, ipoint )
      {
        if ( gamma[ipoint] > 0 )
        {
          motion_[motionCounter] = mpatch[ipoint];
          motionCounter++;
        }
      }
    }
  }

  void RBFMotionSolverExt::setMotion( const vectorField & m )
  {
    if ( m.size() != movingIDs_.size() )
    {
      FatalErrorIn( "void RBFMotionSolverExt::setMotion(const vectorField& m)" )
      << "Incorrect size of motion points: m = " << m.size()
      << " movingIDs = " << movingIDs_.size()
      << abort( FatalError );
    }

    motion_ = vector::zero;

    forAll( m, i )
    {
      motion_[i] = m[i];
    }
  }

  const vectorField & RBFMotionSolverExt::movingPoints() const
  {
    return movingPoints_;
  }

  tmp<pointField> RBFMotionSolverExt::curPoints() const
  {
    // Prepare new points: same as old point
    tmp<pointField> tnewPoints
    (
      new vectorField( mesh().nPoints(), vector::zero )
    );

    pointField & newPoints = tnewPoints();

    // Add motion to newPoints
    forAll( motion_, i )
    {
      newPoints[movingIDs_[i]] = motion_[i];
    }

    // Write the new point discplacement to newPoints. In newPoints movement of all (local) moving points is already included
    moveMesh( newPoints );

    // Set the other side of the TwoD mesh is required
    if ( twoDPointCorrectorRBF_.required() )
    {
      twoDPointCorrectorRBF_.setShadowSide( newPoints );
    }

    // 4. Add old point positions
    newPoints += mesh().points();

    return tnewPoints;
  }

  void RBFMotionSolverExt::solve()
  {
    // Motion must be set by using setMotion();
  }

  void RBFMotionSolverExt::updateMesh( const mapPolyMesh & )
  {
    // Recalculate control point IDs
    makeControlIDs();
  }
}

// ************************************************************************* //
