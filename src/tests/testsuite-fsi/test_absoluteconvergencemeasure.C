
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "AbsoluteConvergenceMeasure.H"
#include "gtest/gtest.h"

using namespace fsi;

TEST(AbsoluteConvergenceMeasure, measure)
{
    scalar convergenceLimit = 9.0;

    AbsoluteConvergenceMeasure measure(0, false, convergenceLimit);

    // Create data sets for old state of data and new state of data
    fsi::vector oldValues0(3), oldValues1(3), oldValues2(3), newValues(3);
    oldValues0 << -2.0, -1.0, 0.0;
    oldValues1 << 2.0, 3.0, 4.0;
    oldValues2 << 3.0, 4.0, 5.0;
    newValues << 5.0, 6.0, 7.0;

    measure.measure(oldValues0, newValues);
    ASSERT_FALSE(measure.isConvergence());

    measure.printState();

    measure.measure(oldValues1, newValues);
    ASSERT_TRUE(measure.isConvergence());

    measure.printState();

    measure.measure(oldValues2, newValues);
    ASSERT_TRUE(measure.isConvergence());

    measure.printState();
}
