
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "MultiLevelFsiSolver.H"

using namespace fsi;

MultiLevelFsiSolver::MultiLevelFsiSolver(
  shared_ptr<MultiLevelSolver> fluidSolver,
  shared_ptr<MultiLevelSolver> solidSolver,
  shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures,
  bool parallel,
  int extrapolationOrder
  )
  :
  FsiSolver( fluidSolver->solver, solidSolver->solver, convergenceMeasures, parallel, extrapolationOrder ),
  fluidSolver( fluidSolver ),
  solidSolver( solidSolver ),
  useJacobian( false )
{
  assert( this->fluid );
  assert( this->solid );
  assert( fluidSolver );
  assert( solidSolver );
  assert( this->convergenceMeasures );
  assert( fluidSolver->couplingGridSize > 0 );
  assert( convergenceMeasures->size() > 0 );
  assert( fluidSolver->level == solidSolver->level );
  assert( fluidSolver->participantId == 0 );
  assert( solidSolver->participantId == 1 );

  fluidSolver->updateInterpolationMatrices();
  solidSolver->updateInterpolationMatrices();

  // Interpolate initial solution from own mesh to coupling grid.

  matrix solidtmp, fluidtmp;

  solidSolver->rbfInterpToCouplingMesh->interpolate( solid->data, solidtmp );
  fluidSolver->rbfInterpToCouplingMesh->interpolate( fluid->data, fluidtmp );

  assert( solidtmp.rows() == solidSolver->couplingGridSize );
  assert( solidtmp.cols() == solid->dim );
  assert( fluidtmp.rows() == fluidSolver->couplingGridSize );
  assert( fluidtmp.cols() == fluid->dim );

  // Initialize control variable x
  if ( parallel )
  {
    x.resize( solidtmp.rows() * solidtmp.cols() + fluidtmp.rows() * fluidtmp.cols() );

    x.head( solidtmp.rows() * solidtmp.cols() ) = Eigen::Map<fsi::vector>( solidtmp.data(), solidtmp.rows() * solidtmp.cols() );
    x.tail( fluidtmp.rows() * fluidtmp.cols() ) = Eigen::Map<fsi::vector>( fluidtmp.data(), fluidtmp.rows() * fluidtmp.cols() );
  }

  if ( !parallel )
    x = Eigen::Map<fsi::vector>( solidtmp.data(), solidtmp.rows() * solidtmp.cols() );

  // Save initial solution for extrapolation
  previousSolutions.clear();
  previousSolutions.push_front( x );

  assert( previousSolutions.size() == 1 );
}

MultiLevelFsiSolver::~MultiLevelFsiSolver()
{}

void MultiLevelFsiSolver::evaluate(
  const fsi::vector & input,
  fsi::vector & output,
  fsi::vector & R
  )
{
  assert( init );
  assert( fluidSolver->solver->init );
  assert( solidSolver->solver->init );
  assert( fluidSolver->solver->t > 0 );
  assert( input.rows() == R.rows() );
  assert( input.rows() == solidSolver->couplingGridSize * solid->dim || input.rows() == solidSolver->couplingGridSize * solid->dim + fluidSolver->couplingGridSize * fluid->dim );
  assert( output.rows() == input.rows() );

  Info << endl << "Time = " << fluidSolver->solver->t << ", iteration = " << iter + 1 << ", level = " << fluidSolver->level << endl;

  matrix a = Eigen::Map<const matrix> ( input.head( solidSolver->couplingGridSize * solid->dim ).data(), solidSolver->couplingGridSize, solid->dim );

  if ( !parallel )
  {
    assert( input.rows() == solidSolver->couplingGridSize * solid->dim );

    // Initialize variables
    matrix p( fluidSolver->couplingGridSize, fluid->dim ), aout( solidSolver->couplingGridSize, solid->dim );

    fluidSolver->solve( a, p );
    solidSolver->solve( p, aout );

    output = Eigen::Map<fsi::vector> ( aout.data(), aout.rows() * aout.cols() );

    fluid->couplingData.setDataOld( fluid->couplingData.data );
    fluid->couplingData.setData( p );
    solid->couplingData.setDataOld( a );
    solid->couplingData.setData( aout );
  }

  if ( parallel )
  {
    assert( input.rows() == solidSolver->couplingGridSize * solid->dim + fluidSolver->couplingGridSize * fluid->dim );

    // Initialize variables
    matrix p( fluidSolver->couplingGridSize, fluid->dim ), pout( fluidSolver->couplingGridSize, fluid->dim ), aout( solidSolver->couplingGridSize, solid->dim );

    p = Eigen::Map<const matrix> ( input.tail( fluidSolver->couplingGridSize * fluid->dim ).data(), fluidSolver->couplingGridSize, fluid->dim );

    fluidSolver->solve( a, pout );
    solidSolver->solve( p, aout );

    output.head( solidSolver->couplingGridSize * solid->dim ) = Eigen::Map<fsi::vector> ( aout.data(), aout.rows() * aout.cols() );
    output.tail( fluidSolver->couplingGridSize * fluid->dim ) = Eigen::Map<fsi::vector> ( pout.data(), pout.rows() * pout.cols() );

    fluid->couplingData.setDataOld( p );
    fluid->couplingData.setData( pout );
    solid->couplingData.setDataOld( a );
    solid->couplingData.setData( aout );
  }

  // Calculate residual
  R = output - input;

  // Increment iterators
  iter++;
  nbIter++;

  x = input;

  if ( useJacobian )
  {
    assert( x.rows() == xf.rows() );
    assert( x.rows() == J.rows() );
    assert( J.rows() == J.cols() );

    Info << "Output space mapping: include Jacobian information in the residual function" << endl;

    R += J * (x - xf);
  }
}

void MultiLevelFsiSolver::initTimeStep()
{
  assert( !init );

  FsiSolver::initTimeStep();
}

bool MultiLevelFsiSolver::isConvergence()
{
  assert( init );
  assert( convergenceMeasures->size() > 0 );
  assert( fluid->t > 0 );

  allConverged = true;

  for ( std::list<std::shared_ptr<ConvergenceMeasure> >::const_iterator iterator = convergenceMeasures->begin(), end = convergenceMeasures->end(); iterator != end; ++iterator )
  {
    // Initialize variables
    std::shared_ptr<ConvergenceMeasure> measure = *iterator;
    fsi::vector oldValues, newValues;

    assert( measure->dataId == 0 || measure->dataId == 1 );

    // Subtract initial values
    if ( measure->dataId == 0 )
    {
      oldValues = Eigen::Map<fsi::vector> ( solid->couplingData.dataprev.data(), solid->couplingData.dataprev.rows() * solid->couplingData.dataprev.cols() );
      newValues = Eigen::Map<fsi::vector> ( solid->couplingData.data.data(), solid->couplingData.data.rows() * solid->couplingData.dataprev.cols() );

      oldValues.array() -= solid->initialValue;
      newValues.array() -= solid->initialValue;
    }

    if ( measure->dataId == 1 )
    {
      oldValues = Eigen::Map<fsi::vector> ( fluid->couplingData.dataprev.data(), fluid->couplingData.dataprev.rows() * fluid->couplingData.dataprev.cols() );
      newValues = Eigen::Map<fsi::vector> ( fluid->couplingData.data.data(), fluid->couplingData.data.rows() * fluid->couplingData.dataprev.cols() );

      oldValues.array() -= fluid->initialValue;
      newValues.array() -= fluid->initialValue;
    }

    // Measure convergence
    measure->measure( oldValues, newValues );

    // Print state
    measure->printState();

    if ( !measure->isConvergence() )
      allConverged = false;
  }

  if ( allConverged )
  {
    Info << "All converged" << endl;
    fluid->couplingData.dataprev.setZero();
    solid->couplingData.dataprev.setZero();
  }

  return allConverged;
}

bool MultiLevelFsiSolver::isConvergence(
  const fsi::vector & x,
  const fsi::vector & xprev
  )
{
  assert( init );
  assert( convergenceMeasures->size() > 0 );
  assert( x.rows() == xprev.rows() );
  assert( x.rows() == solidSolver->couplingGridSize * solid->dim || x.rows() == solidSolver->couplingGridSize * solid->dim + fluidSolver->couplingGridSize * fluid->dim );
  assert( fluid->t > 0 );

  allConverged = true;

  for ( std::list<std::shared_ptr<ConvergenceMeasure> >::const_iterator iterator = convergenceMeasures->begin(), end = convergenceMeasures->end(); iterator != end; ++iterator )
  {
    // Initialize variables
    std::shared_ptr<ConvergenceMeasure> measure = *iterator;
    fsi::vector oldValues, newValues;

    // Subtract initial values
    if ( measure->dataId == 0 )
    {
      oldValues = xprev.head( solidSolver->couplingGridSize * solid->dim ).array() - solid->initialValue;
      newValues = x.head( solidSolver->couplingGridSize * solid->dim ).array() - solid->initialValue;
    }

    if ( measure->dataId == 1 )
    {
      if ( parallel )
      {
        assert( xprev.rows() == fluidSolver->couplingGridSize * fluid->dim + solidSolver->couplingGridSize * solid->dim );

        oldValues = xprev.tail( fluidSolver->couplingGridSize * fluid->dim );
        newValues = x.tail( fluidSolver->couplingGridSize * fluid->dim );
      }

      if ( !parallel )
      {
        oldValues = Eigen::Map<fsi::vector> ( fluid->couplingData.dataprev.data(), fluid->couplingData.dataprev.rows() * fluid->couplingData.dataprev.cols() );
        newValues = Eigen::Map<fsi::vector> ( fluid->couplingData.data.data(), fluid->couplingData.data.rows() * fluid->couplingData.dataprev.cols() );
      }

      oldValues.array() -= fluid->initialValue;
      newValues.array() -= fluid->initialValue;
    }

    // Measure convergence
    measure->measure( oldValues, newValues );

    // Print state
    measure->printState();

    if ( !measure->isConvergence() )
      allConverged = false;
  }

  if ( allConverged )
  {
    Info << "All converged" << endl;
    fluid->couplingData.dataprev.setZero();
    solid->couplingData.dataprev.setZero();
  }

  return allConverged;
}

void MultiLevelFsiSolver::setSurrogateData(
  fsi::vector & xf,
  matrix & J
  )
{
  this->xf = xf;
  this->J = J;
}

void MultiLevelFsiSolver::setUseJacobian( bool useJacobian )
{
  this->useJacobian = useJacobian;
}
