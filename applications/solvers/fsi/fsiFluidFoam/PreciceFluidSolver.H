
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

  void readDataAcousticsTemperatureGradient();

  void readDataAcousticsPressureGradient();

  void readDataAcousticsVelocityGradient();

  void setReadPositions();

  void setReadPositionsAcoustics();

  void setWritePositions();

  void setWritePositionsAcoustics();

  void writeData( const matrix & data );

  void writeDataAcoustics();

  void writeAcousticsDensity();

  void writeAcousticsPressure();

  void writeAcousticsVelocity();

  shared_ptr<foamFluidSolver> solver;
  shared_ptr<precice::SolverInterface> precice;

  Eigen::VectorXi idsReadPositions;
  Eigen::VectorXi idsWritePositions;
  Eigen::VectorXi idsWritePositionsAcoustics;

  const std::string FLUID_ACOUSTICS_WRITE;
  const std::string FLUID_ACOUSTICS_READ;
};

#endif
