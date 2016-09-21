
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef PreciceSolidSolver_H
#define PreciceSolidSolver_H

#include <precice/SolverInterface.hpp>

#include "foamSolidSolver.H"
#include "Solver.H"

using std::shared_ptr;

using namespace fsi;

class PreciceSolidSolver : public Solver
{
public:

  explicit PreciceSolidSolver( shared_ptr<foamSolidSolver> solver );

  virtual ~PreciceSolidSolver();

  virtual void run();

  void init();

  void readData( matrix & data );

  void setReadPositions();

  void setWritePositions();

  void writeData( const matrix & data );

  shared_ptr<foamSolidSolver> solver;
  shared_ptr<precice::SolverInterface> precice;

  Eigen::VectorXi idsReadPositions;
  Eigen::VectorXi idsWritePositions;
};

#endif
