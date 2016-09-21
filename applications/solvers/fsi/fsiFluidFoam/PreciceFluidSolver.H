
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef PreciceFluidSolver_H
#define PreciceFluidSolver_H

#include <precice/SolverInterface.hpp>

#include "foamFluidSolver.H"
#include "Solver.H"

using namespace fsi;

class PreciceFluidSolver : public Solver
{
public:

  explicit PreciceFluidSolver( shared_ptr<foamFluidSolver> solver );

  virtual ~PreciceFluidSolver();

  virtual void run();

  void init();

  void readData( matrix & data );

  void readDataAcoustics();

  void setReadPositions();

  void setWritePositions();

  void setWritePositionsAcoustics();

  void writeData( const matrix & data );

  void writeDataAcoustics();

  void writeAcousticsDensity();

  void writeAcousticsPressure();

  void writeAcousticsVelocity();

  void writeAcousticsDensityGradient();

  void writeAcousticsPressureGradient();

  void writeAcousticsVelocityGradient();

  shared_ptr<foamFluidSolver> solver;
  shared_ptr<precice::SolverInterface> precice;

  Eigen::VectorXi idsReadPositions;
  Eigen::VectorXi idsWritePositions;
  Eigen::VectorXi idsWritePositionsAcoustics;
};

#endif
