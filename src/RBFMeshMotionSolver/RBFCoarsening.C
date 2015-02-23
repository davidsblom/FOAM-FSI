
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
    bool livePointSelection,
    double tol,
    double tolLivePointSelection,
    int coarseningMinPoints,
    int coarseningMaxPoints
    )
    :
    rbf( rbf ),
    rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction ) ) ),
    enabled( enabled ),
    livePointSelection( livePointSelection ),
    tol( tol ),
    tolLivePointSelection( tolLivePointSelection ),
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

  RBFCoarsening::RBFCoarsening(
    std::shared_ptr<RBFInterpolation> rbf,
    bool enabled,
    bool livePointSelection,
    double tol,
    int coarseningMinPoints,
    int coarseningMaxPoints
    )
    :
    rbf( rbf ),
    rbfCoarse( std::shared_ptr<RBFInterpolation> ( new RBFInterpolation( rbf->rbfFunction ) ) ),
    enabled( enabled ),
    livePointSelection( livePointSelection ),
    tol( tol ),
    tolLivePointSelection( tol ),
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
    this->positions = positions;
    this->positionsInterpolation = positionsInterpolation;
  }

  void RBFCoarsening::greedySelection( const matrix & values )
  {
    assert( positions.cols() == positionsInterpolation.cols() );
    assert( positions.cols() > 0 );
    assert( positions.rows() > 0 );
    assert( positionsInterpolation.rows() > 0 );
    assert( positions.rows() == values.rows() );

    matrix usedPositions = positions;

    if ( enabled )
    {
      // Greedy algorithm

      rbf::vector errorList( positions.rows() );
      selectedPositions.resize( 1 );
      selectedPositions.setZero();
      rbf::matrix positionsInterpolationCoarse = positions;

      int maxNbPoints = std::min( coarseningMaxPoints, static_cast<int>( positions.rows() ) );
      int minPoints = std::min( coarseningMinPoints, static_cast<int>( positions.rows() ) );
      double error = 0;

      // Create RBF interpolator
      rbfCoarse->Phi.resize( 0, 0 );

      for ( int i = selectedPositions.rows(); i < maxNbPoints + 1; i++ )
      {
        // Build the matrices used for the RBF interpolation
        rbf::matrix positionsCoarse( i, positions.cols() );
        rbf::matrix valuesCoarse( positionsCoarse.rows(), positionsCoarse.cols() );
        rbf::matrix valuesInterpolationCoarse( positionsInterpolationCoarse.rows(), positionsInterpolationCoarse.cols() );

        for ( int j = 0; j < selectedPositions.rows(); j++ )
        {
          positionsCoarse.row( j ) = positions.row( selectedPositions( j ) );
          valuesCoarse.row( j ) = values.row( selectedPositions( j ) );
        }

        // Perform the RBF interpolation.
        rbfCoarse->interpolate( positionsCoarse, positionsInterpolationCoarse, valuesCoarse, valuesInterpolationCoarse );

        // Evaluate the error
        for ( int j = 0; j < valuesInterpolationCoarse.rows(); j++ )
          errorList( j ) = ( valuesInterpolationCoarse.row( j ).array() - values.row( j ).array() ).matrix().norm();

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

        double epsilon = std::sqrt( SMALL );
        error = ( valuesInterpolationCoarse.array() - values.array() ).matrix().norm() / (values.norm() + epsilon);
        bool convergence = (error < tol && i >= minPoints) || i >= maxNbPoints;

        if ( convergence )
          break;

        selectedPositions.conservativeResize( selectedPositions.rows() + 1 );
        selectedPositions( selectedPositions.rows() - 1 ) = index;
      }

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
      if ( livePointSelection )
      {
        // For RBF mesh interpolation, the values to be interpolated need to be
        // the total displacements. As input, the incremental displacements
        // are given.
        if ( this->values.cols() != values.cols() )
          this->values = values;
        else
          this->values.array() += values.array();

        // Check if re-selection is necessary
        bool reselection = true;

        if ( rbfCoarse->computed )
        {
          rbf::matrix valuesCoarse( selectedPositions.rows(), values.cols() );
          rbf::matrix valuesInterpolationCoarse( positions.rows(), valuesInterpolation.cols() );

          for ( int j = 0; j < selectedPositions.rows(); j++ )
            valuesCoarse.row( j ) = this->values.row( selectedPositions( j ) );

          rbfCoarse->interpolate2( valuesCoarse, valuesInterpolationCoarse );

          double epsilon = std::sqrt( SMALL );
          double error = ( valuesInterpolationCoarse.array() - this->values.array() ).matrix().norm() / (this->values.norm() + epsilon);
          bool convergence = error < tolLivePointSelection;

          if ( convergence )
            reselection = false;

          Info << "RBF interpolation coarsening: error = " << error << ", tol = " << tolLivePointSelection << ", reselection = ";

          if ( reselection )
            Info << "true";
          else
            Info << "false";

          Info << endl;
        }

        if ( reselection )
        {
          greedySelection( this->values );

          rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
        }
      }
      else
      if ( !rbf->computed )
      {
        // Unit displacement of control points
        matrix unitDisplacement( positions.rows(), positions.cols() );
        unitDisplacement.fill( 1 );
        greedySelection( unitDisplacement );

        rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
      }

      rbf::matrix selectedValues( selectedPositions.rows(), values.cols() );

      for ( int i = 0; i < selectedValues.rows(); i++ )
        selectedValues.row( i ) = values.row( selectedPositions( i ) );

      usedValues = selectedValues;
    }
    else
    {
      if ( !rbf->computed )
      {
        rbf->compute( positions, positionsInterpolation );
        rbf->Hhat.conservativeResize( rbf->Hhat.rows(), rbf->Hhat.cols() - nbStaticFaceCentersRemove );
      }
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
