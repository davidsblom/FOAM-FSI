
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFCoarsening.H"

namespace rbf
{
  RBFCoarsening::RBFCoarsening(
    std::shared_ptr<RBFInterpolation> rbf,
    bool enabled,
    double tol,
    int coarseningMinPoints,
    int coarseningMaxPoints
    )
    :
    rbf( rbf ),
    enabled( enabled ),
    tol( tol ),
    coarseningMinPoints( coarseningMinPoints ),
    coarseningMaxPoints( coarseningMaxPoints ),
    selectedPositions( 1 ),
    nbStaticFaceCentersRemove( 0 )
  {
    assert( rbf );
    assert( coarseningMinPoints <= coarseningMaxPoints );
    assert( coarseningMinPoints > 0 );
    assert( coarseningMaxPoints > 0 );
    assert( tol > 0 );
    assert( tol < 1 );
  }

  /* Select a subset of control point with a greedy algorithm.
   * The selection of the points is based on a displacement/motion of
   * a unit displacement at every control point. Based on a user specified tolerance,
   * the control points are selected.
   */
  void RBFCoarsening::compute(
    const matrix & positions,
    const matrix & positionsInterpolation
    )
  {
    matrix usedPositions = positions;

    if ( enabled )
    {
      // Greedy algorithm

      rbf::vector errorList( positions.rows() );
      selectedPositions.setZero();
      rbf::matrix positionsInterpolationCoarse = positions;

      int maxNbPoints = std::min( coarseningMaxPoints, static_cast<int>( positions.rows() ) );
      int minPoints = std::min( coarseningMinPoints, static_cast<int>( positions.rows() ) );
      double error = 0;

      // Create RBF interpolator
      rbf::RBFInterpolation rbfCoarse( rbf->rbfFunction );

      for ( int i = 1; i < maxNbPoints; i++ )
      {
        // Build the matrices used for the RBF interpolation
        rbf::matrix positionsCoarse( i, positions.cols() );

        for ( int j = 0; j < selectedPositions.rows(); j++ )
          positionsCoarse.row( j ) = positions.row( selectedPositions( j ) );

        rbf::matrix valuesCoarse( positionsCoarse.rows(), positionsCoarse.cols() ), valuesInterpolationCoarse( positionsInterpolationCoarse.rows(), positionsInterpolationCoarse.cols() );
        valuesCoarse.setZero();

        valuesCoarse.fill( 1 );

        // Perform the RBF interpolation.
        rbfCoarse.interpolate( positionsCoarse, positionsInterpolationCoarse, valuesCoarse, valuesInterpolationCoarse );

        // Evaluate the error
        for ( int j = 0; j < valuesInterpolationCoarse.rows(); j++ )
          errorList( j ) = (valuesInterpolationCoarse.row( j ).array() - 1).matrix().norm();

        // Select the point with the largest error which is not already selected.

        int index = -1;
        double largestError = -1;

        for ( int j = 0; j < errorList.rows(); j++ )
        {
          bool notSelected = true;

          for ( int k = 0; k < selectedPositions.rows(); k++ )
            if ( selectedPositions( k ) == j )
              notSelected = false;

          if ( errorList( j ) > largestError && notSelected )
          {
            index = j;
            largestError = errorList( j );
          }
        }

        double N = valuesInterpolationCoarse.rows() * valuesInterpolationCoarse.cols();
        error = std::sqrt( (valuesInterpolationCoarse.array() - 1).matrix().squaredNorm() / N );
        bool convergence = error < tol && i >= minPoints;

        if ( convergence )
          break;

        selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
        selectedPositions( selectedPositions.rows() - 1 ) = index;
      }

      std::sort( selectedPositions.data(), selectedPositions.data() + selectedPositions.size() );

      Info << "RBF interpolation coarsening: selected " << selectedPositions.rows() << "/" << positions.rows() << " points, error = " << error << ", tol = " << tol << endl;

      rbf::matrix positionsCoarse( selectedPositions.rows(), positions.cols() );

      for ( int i = 0; i < selectedPositions.rows(); i++ )
        positionsCoarse.row( i ) = positions.row( selectedPositions( i ) );

      usedPositions = positionsCoarse;
    }

    rbf->compute( usedPositions, positionsInterpolation );
  }

  void RBFCoarsening::computeInParallel(
    const matrix & positions,
    const matrix & positionsInterpolation,
    const Eigen::VectorXi & positionsParallelLocation
    )
  {
    assert( positions.rows() == positionsParallelLocation.rows() );

    matrix usedPositions = positions;

    if ( enabled )
    {
      // Parallel Greedy algorithm

      selectedPositions.setZero();

      // Build the positionsInterpolationCoarse matrix based on positionsParallelLocation

      int nbLocalPoints = 0;

      for ( int i = 0; i < positions.rows(); i++ )
        if ( positionsParallelLocation( i ) == 1 )
          nbLocalPoints++;

      assert( nbLocalPoints <= positions.rows() );

      // Build the reverse table of positionsParallelLocation
      Eigen::VectorXi positionsParallelLocationReverse( nbLocalPoints );
      rbf::vector errorList( nbLocalPoints );
      int index = 0;

      for ( int i = 0; i < positions.rows(); i++ )
      {
        if ( positionsParallelLocation( i ) == 1 )
        {
          positionsParallelLocationReverse( index ) = i;
          index++;
        }
      }

      assert( index == nbLocalPoints );

      // Build the local positions matrix
      rbf::matrix positionsInterpolationCoarse( nbLocalPoints, positions.cols() );

      for ( int i = 0; i < positionsInterpolationCoarse.rows(); i++ )
        positionsInterpolationCoarse.row( i ) = positions.row( positionsParallelLocationReverse( i ) );

      // Start the Greedy algorithm based on the local positions.

      int maxNbPoints = std::min( coarseningMaxPoints, static_cast<int>( positions.rows() ) );
      int minPoints = std::min( coarseningMinPoints, static_cast<int>( positions.rows() ) );
      double error = 0;

      // Create RBF interpolator
      rbf::RBFInterpolation rbfCoarse( rbf->rbfFunction );

      for ( int i = 1; i < maxNbPoints; i++ )
      {
        // Build the matrices used for the RBF interpolation
        rbf::matrix positionsCoarse( i, positions.cols() );

        for ( int j = 0; j < selectedPositions.rows(); j++ )
          positionsCoarse.row( j ) = positions.row( selectedPositions( j ) );

        rbf::matrix valuesCoarse( positionsCoarse.rows(), positionsCoarse.cols() ), valuesInterpolationCoarse( positionsInterpolationCoarse.rows(), positionsInterpolationCoarse.cols() );
        valuesCoarse.setZero();

        valuesCoarse.fill( 1 );

        // Perform the RBF interpolation.
        rbfCoarse.interpolate( positionsCoarse, positionsInterpolationCoarse, valuesCoarse, valuesInterpolationCoarse );

        // Evaluate the error
        for ( int j = 0; j < valuesInterpolationCoarse.rows(); j++ )
          errorList( j ) = (valuesInterpolationCoarse.row( j ).array() - 1).matrix().norm();

        // Select the point with the largest error which is not already selected.

        int index = -1;
        double largestError = -1;

        for ( int j = 0; j < errorList.rows(); j++ )
        {
          bool notSelected = true;

          for ( int k = 0; k < selectedPositions.rows(); k++ )
            if ( selectedPositions( k ) == j )
              notSelected = false;

          if ( errorList( j ) > largestError && notSelected )
          {
            index = j;
            largestError = errorList( j );
          }
        }

        // Communicate the index in the positions matrix of the largest error,
        // and the error to every processor.
        labelList globalErrorIndices( Pstream::nProcs(), 0 );
        scalarList globalErrorList( Pstream::nProcs(), 0.0 );
        globalErrorIndices[Pstream::myProcNo()] = positionsParallelLocationReverse( index );
        globalErrorList[Pstream::myProcNo()] = largestError;

        reduce( globalErrorIndices, sumOp<labelList>() );
        reduce( globalErrorList, sumOp<scalarList>() );

        // Select the positions with the largest of the global surface mesh.
        // The point with the largest error is automatically not selected
        // already.
        index = -1;
        largestError = -1;

        for ( int j = 0; j < globalErrorList.size(); j++ )
        {
          if ( globalErrorList[j] > largestError )
          {
            index = globalErrorIndices[j];
            largestError = globalErrorList[j];
          }
        }

        // absolute error < tol.
        // This is justified for an applied unit displacement.
        error = largestError;
        bool convergence = error < tol && i >= minPoints;

        if ( convergence )
          break;

        selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
        selectedPositions( selectedPositions.rows() - 1 ) = index;
      }

      std::sort( selectedPositions.data(), selectedPositions.data() + selectedPositions.size() );

      Info << "RBF interpolation coarsening: selected " << selectedPositions.rows() << "/" << positions.rows() << " points, error = " << error << ", tol = " << tol << endl;

      rbf::matrix positionsCoarse( selectedPositions.rows(), positions.cols() );

      for ( int i = 0; i < selectedPositions.rows(); i++ )
        positionsCoarse.row( i ) = positions.row( selectedPositions( i ) );

      usedPositions = positionsCoarse;
    }

    rbf->compute( usedPositions, positionsInterpolation );
  }

  void RBFCoarsening::interpolate(
    const matrix & values,
    matrix & valuesInterpolation
    )
  {
    matrix usedValues = values;

    if ( enabled )
    {
      rbf::matrix selectedValues( selectedPositions.rows(), values.cols() );

      for ( int i = 0; i < selectedValues.rows(); i++ )
        selectedValues.row( i ) = values.row( selectedPositions( i ) );

      usedValues = selectedValues;
    }

    usedValues.conservativeResize( usedValues.rows() - nbStaticFaceCentersRemove, usedValues.cols() );

    rbf->interpolate( usedValues, valuesInterpolation );
  }

  void RBFCoarsening::setNbMovingAndStaticFaceCenters(
    int nbMovingFaceCenters,
    int nbStaticFaceCenters
    )
  {
    nbStaticFaceCentersRemove = nbStaticFaceCenters;

    if ( enabled )
    {
      // Determine the number of selected static face centers
      nbStaticFaceCentersRemove = 0;

      for ( int i = 0; i < selectedPositions.rows(); i++ )
        if ( selectedPositions( i ) >= nbMovingFaceCenters )
          nbStaticFaceCentersRemove++;
    }

    rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
  }
}
