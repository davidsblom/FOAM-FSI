
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "fvCFD.H"
#include "MinIterationConvergenceMeasure.H"

namespace fsi
{
  MinIterationConvergenceMeasure::MinIterationConvergenceMeasure(
    int dataId,
    int minimumIterationCount
    )
    :
    ConvergenceMeasure( dataId, minimumIterationCount ),
    minimumIterationCount( minimumIterationCount ),
    currentIteration( 0 )
  {
    assert( minimumIterationCount >= 0 );
  }

  void MinIterationConvergenceMeasure::measure(
    vector & oldValues,
    vector & newValues
    )
  {
    assert( oldValues.rows() == newValues.rows() );

    currentIteration++;

    isConvergence_ = minimumIterationCount <= currentIteration;
  }

  void MinIterationConvergenceMeasure::newMeasurementSeries()
  {
    currentIteration = 0;
    isConvergence_ = false;
  }

  bool MinIterationConvergenceMeasure::isConvergence()
  {
    return isConvergence_;
  }

  void MinIterationConvergenceMeasure::printState()
  {
    Info << "min iteration convergence measure: #it = "
         << currentIteration
         << " of "
         << minimumIterationCount
         << ", conv = ";

    if ( isConvergence_ )
      Info << "true";

    if ( !isConvergence_ )
      Info << "false";

    Info << endl;
  }
}
