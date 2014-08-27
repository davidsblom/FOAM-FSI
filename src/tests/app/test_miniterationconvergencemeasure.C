
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MinIterationConvergenceMeasure.H"
#include "gtest/gtest.h"

using namespace fsi;

TEST( MinIterationConvergenceMeasure, measure )
{
  int minimumIterationCount = 3;

  MinIterationConvergenceMeasure measure( 0, minimumIterationCount );

  // Create data sets for old state of data and new state of data
  vector oldValues0( 3 ), oldValues1( 3 ), oldValues2( 3 ), newValues( 3 );
  oldValues0 << 1, 1, 1;
  oldValues1 << 2, 2, 2;
  oldValues2 << 2.9, 2.9, 2.9;
  newValues << 3, 3, 3;

  measure.measure( oldValues0, newValues );
  ASSERT_FALSE( measure.isConvergence() );

  measure.measure( oldValues1, newValues );
  ASSERT_FALSE( measure.isConvergence() );

  measure.measure( oldValues2, newValues );
  ASSERT_TRUE( measure.isConvergence() );

  measure.measure( oldValues2, newValues );
  ASSERT_TRUE( measure.isConvergence() );

  measure.printState();
}

TEST( MinIterationConvergenceMeasure, newMeasurementSeries )
{
  int minimumIterationCount = 3;

  MinIterationConvergenceMeasure measure( 0, minimumIterationCount );

  // Create data sets for old state of data and new state of data
  vector oldValues0( 3 ), oldValues1( 3 ), oldValues2( 3 ), newValues( 3 );
  oldValues0 << 1, 1, 1;
  oldValues1 << 2, 2, 2;
  oldValues2 << 2.9, 2.9, 2.9;
  newValues << 3, 3, 3;

  measure.measure( oldValues0, newValues );
  ASSERT_FALSE( measure.isConvergence() );

  measure.measure( oldValues1, newValues );
  ASSERT_FALSE( measure.isConvergence() );

  ASSERT_EQ( measure.currentIteration, 2 );

  measure.newMeasurementSeries();

  ASSERT_EQ( measure.currentIteration, 0 );
}
