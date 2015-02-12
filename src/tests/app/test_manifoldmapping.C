
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "ManifoldMapping.H"
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

using namespace fsi;

class FineModelParabola : public SurrogateModel
{
public:

  FineModelParabola ( double tol )
    :
    tol( tol ),
    iter( 0 ),
    allConverged_( false )
  {
    assert( tol > 0 );
    assert( tol < 1 );
  }

  virtual bool allConverged()
  {
    return allConverged_;
  }

  virtual void evaluate(
    const fsi::vector & x,
    fsi::vector & output,
    fsi::vector & R
    )
  {
    assert( x.rows() == 2 );
    assert( R.rows() == 3 );

    fsi::vector t( 3 );
    t << -1, 0, 1;

    R = x( 0 ) * (x( 1 ) * t.array() + 1.0).pow( 2 );
    output = R;
  }

  virtual bool isConvergence()
  {
    return true;
  }

  virtual bool isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    )
  {
    iter++;

    allConverged_ = false;

    if ( (x - xprev).norm() < tol * x.norm() )
      allConverged_ = true;

    return allConverged_;
  }

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk
    );

  virtual void optimize(
    const fsi::vector & x0,
    fsi::vector & xk
    )
  {}

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk,
    const matrix & B,
    const fsi::vector & xktilde,
    const fsi::vector & xkp
    )
  {}

  virtual void setSurrogateData(
    fsi::vector & xf,
    matrix & J
    )
  {}

  virtual void setUseJacobian( bool useJacobian )
  {}

  double tol;
  int iter;
  bool allConverged_;
};

// Generic functor
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
  typedef _Scalar Scalar;
  enum {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
  typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

  int m_inputs, m_values;

  Functor() : m_inputs( InputsAtCompileTime ), m_values( ValuesAtCompileTime ) {}

  Functor(
    int inputs,
    int values
    ) : m_inputs( inputs ), m_values( values ) {}

  int inputs() const
  {
    return m_inputs;
  }

  int values() const
  {
    return m_values;
  }
};

struct fineModelParabolaFunctor : Functor<double>
{
  fineModelParabolaFunctor(
    FineModelParabola * model,
    const fsi::vector * y
    )
    :
    Functor<double>( 2, 3 ),
    model( model ),
    y( y )
  {}

  int operator()(
    const Eigen::VectorXd & x,
    Eigen::VectorXd & fvec
    ) const
  {
    fsi::vector output( x.rows() );
    model->evaluate( x, output, fvec );
    fvec -= *y;
    return 0;
  }

  FineModelParabola * model;
  const fsi::vector * y;
};

void FineModelParabola::optimize(
  const fsi::vector & y,
  const fsi::vector & x0,
  fsi::vector & xk
  )
{
  fineModelParabolaFunctor functor( this, &y );
  Eigen::NumericalDiff<fineModelParabolaFunctor, Eigen::Central> numDiff( functor );
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<fineModelParabolaFunctor, Eigen::Central>, double> lm( numDiff );

  lm.parameters.maxfev = 2000;
  lm.parameters.xtol = 1.0e-13;
  lm.parameters.ftol = 1.0e-13;

  xk = x0;

  int ret = lm.minimize( xk );

  // 2: RelativeErrorTooSmall
  assert( ret == 2 );
}

class FineModelSimple : public SurrogateModel
{
public:

  FineModelSimple( double tol )
    :
    tol( tol ),
    iter( 0 ),
    allConverged_( false )
  {
    assert( tol > 0 );
    assert( tol < 1 );
  }

  virtual bool allConverged()
  {
    return allConverged_;
  }

  virtual void evaluate(
    const fsi::vector & x,
    fsi::vector & output,
    fsi::vector & R
    )
  {
    assert( x.rows() == 2 );
    assert( R.rows() == 3 );

    fsi::vector t( 3 );
    t << -1, 0, 1;

    R = x( 0 ) * x( 0 ) * t.array() * t.array() + x( 0 ) * t.array() + x( 1 );
    output = R;
  }

  virtual bool isConvergence()
  {
    return true;
  }

  virtual bool isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    )
  {
    assert( x.rows() == xprev.rows() );

    iter++;

    allConverged_ = false;

    if ( (x - xprev).norm() < tol * x.norm() )
      allConverged_ = true;

    return allConverged_;
  }

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk
    );

  virtual void optimize(
    const fsi::vector & x0,
    fsi::vector & xk
    )
  {}

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk,
    const matrix & B,
    const fsi::vector & xktilde,
    const fsi::vector & xkp
    )
  {}

  virtual void setSurrogateData(
    fsi::vector & xf,
    matrix & J
    )
  {}

  virtual void setUseJacobian( bool useJacobian )
  {}

  double tol;
  int iter;
  bool allConverged_;
};

struct fineModelSimpleFunctor : Functor<double>
{
  fineModelSimpleFunctor(
    FineModelSimple * model,
    const fsi::vector * y
    )
    :
    Functor<double>( 2, 3 ),
    model( model ),
    y( y )
  {}

  int operator()(
    const Eigen::VectorXd & x,
    Eigen::VectorXd & fvec
    ) const
  {
    fsi::vector output( x.rows() );
    model->evaluate( x, output, fvec );
    fvec -= *y;
    return 0;
  }

  FineModelSimple * model;
  const fsi::vector * y;
};

void FineModelSimple::optimize(
  const fsi::vector & y,
  const fsi::vector & x0,
  fsi::vector & xk
  )
{
  fineModelSimpleFunctor functor( this, &y );
  Eigen::NumericalDiff<fineModelSimpleFunctor, Eigen::Central> numDiff( functor );
  Eigen::LevenbergMarquardt<Eigen::NumericalDiff<fineModelSimpleFunctor, Eigen::Central>, double> lm( numDiff );

  lm.parameters.maxfev = 2000;
  lm.parameters.xtol = 1.0e-13;
  lm.parameters.ftol = 1.0e-13;

  xk = x0;

  lm.minimize( xk );
}

class CoarseModelParabola : public SurrogateModel
{
public:

  CoarseModelParabola()
  {}

  virtual bool allConverged()
  {
    return true;
  }

  virtual void evaluate(
    const fsi::vector & x,
    fsi::vector & output,
    fsi::vector & R
    )
  {
    assert( x.rows() == 2 );
    assert( R.rows() == 3 );

    fsi::vector t( 3 );
    t << -1, 0, 1;

    R = x( 0 ) * t + x( 1 ) * Eigen::MatrixXd::Ones( 3, 1 );
    output = R;
  }

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk
    )
  {
    assert( y.rows() == 3 );
    assert( x0.rows() == 2 );
    assert( xk.rows() == 2 );

    // Analytical optimum
    xk( 0 ) = (y[2] - y[0]) / 2.0;
    xk( 1 ) = (y[0] + y[1] + y[2]) / 3.0;
  }

  virtual void optimize(
    const fsi::vector & x0,
    fsi::vector & xk
    )
  {
    fsi::vector y( 3 );
    y.setZero();
    optimize( y, x0, xk );
  }

  virtual bool isConvergence()
  {
    return true;
  }

  virtual bool isConvergence(
    const fsi::vector & x,
    const fsi::vector & xprev
    )
  {
    return true;
  }

  virtual void optimize(
    const fsi::vector & y,
    const fsi::vector & x0,
    fsi::vector & xk,
    const matrix & B,
    const fsi::vector & xktilde,
    const fsi::vector & xkp
    )
  {}

  virtual void setSurrogateData(
    fsi::vector & xf,
    matrix & J
    )
  {}

  virtual void setUseJacobian( bool useJacobian )
  {}
};

TEST( ManifoldMapping, specification1 )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelParabola( 1.0e-10 ) );
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  // Design specification
  fsi::vector y( 3 );
  y << 0.081, 0.1, 0.121;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  // Settings
  int maxIter = 200;
  double singularityLimit = 1.0e-15;
  int nbReuse = 0;
  int reuseInformationStartingFromTimeIndex = 0;
  bool updateJacobian = false;

  // Create manifold mapping object
  ManifoldMapping mm = ManifoldMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian );

  mm.performPostProcessing( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), 0.1, 1.0e-10 );
  ASSERT_NEAR( xk( 1 ), 0.1, 1.0e-10 );
  ASSERT_EQ( mm.iter, 14 );
}

TEST( ManifoldMapping, specification2 )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelParabola( 1.0e-11 ) );
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  // Design specification
  fsi::vector y( 3 );
  y << 0.10011, 0.10125, 0.10241;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  // Settings
  int maxIter = 500;
  double singularityLimit = 1.0e-15;
  int nbReuse = 0;
  int reuseInformationStartingFromTimeIndex = 0;
  bool updateJacobian = false;

  // Create manifold mapping object
  ManifoldMapping mm = ManifoldMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian );

  mm.performPostProcessing( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), 0.101, 1.0e-3 );
  ASSERT_NEAR( xk( 1 ), 0.006, 1.0e-3 );
  ASSERT_EQ( mm.iter, 322 );
}

TEST( ManifoldMapping, specification3 )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelParabola( 1.0e-7 ) );
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  // Design specification
  fsi::vector y( 3 );
  y << 0, -0.4, 0.1;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  // Settings
  int maxIter = 200;
  double singularityLimit = 1.0e-15;
  int nbReuse = 0;
  int reuseInformationStartingFromTimeIndex = 0;
  bool updateJacobian = false;

  // Create manifold mapping object
  ManifoldMapping mm = ManifoldMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian );

  mm.performPostProcessing( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), -0.101, 1.0e-3 );
  ASSERT_NEAR( xk( 1 ), -0.141, 1.0e-3 );
  ASSERT_EQ( mm.iter, 63 );
}

TEST( ManifoldMapping, specification4 )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelParabola( 1.0e-7 ) );
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  // Design specification
  fsi::vector y( 3 );
  y << 0, -0.35, 0.2;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  // Settings
  int maxIter = 2000;
  double singularityLimit = 1.0e-15;
  int nbReuse = 0;
  int reuseInformationStartingFromTimeIndex = 0;
  bool updateJacobian = false;

  // Create manifold mapping object
  ManifoldMapping mm = ManifoldMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian );

  mm.performPostProcessing( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), 0.007, 1.0e-3 );
  ASSERT_NEAR( xk( 1 ), 4.007, 1.0e-3 );
  ASSERT_EQ( mm.iter, 552 );
}

TEST( ManifoldMapping, specification5 )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelSimple( 1.0e-10 ) );
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  // Design specification
  fsi::vector y( 3 );
  y << 0.081, 0.1, 0.121;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  // Settings
  int maxIter = 200;
  double singularityLimit = 1.0e-15;
  int nbReuse = 0;
  int reuseInformationStartingFromTimeIndex = 0;
  bool updateJacobian = false;

  // Create manifold mapping object
  ManifoldMapping mm = ManifoldMapping( fineModel, coarseModel, maxIter, nbReuse, reuseInformationStartingFromTimeIndex, singularityLimit, updateJacobian );

  mm.performPostProcessing( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), 0.02000799, 1.0e-5 );
  ASSERT_NEAR( xk( 1 ), 0.10039978, 1.0e-5 );
  ASSERT_EQ( mm.iter, 3 );
}

TEST( ManifoldMapping, fineModelSimpleOptimum )
{
  shared_ptr<FineModelSimple> fineModel( new FineModelSimple( 1.0e-7 ) );

  // Design specification
  fsi::vector y( 3 );
  y << 0.081, 0.1, 0.121;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  fineModel->optimize( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), 0.02000799, 1.0e-5 );
  ASSERT_NEAR( xk( 1 ), 0.10039978, 1.0e-5 );
}

TEST( ManifoldMapping, fineModelParabolaOptimum )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelParabola( 1.0e-7 ) );

  // Design specification
  fsi::vector y( 3 );
  y << 0.081, 0.1, 0.121;

  // Initial solution
  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();
  xk.setZero();

  fineModel->optimize( y, x0, xk );

  ASSERT_EQ( xk.rows(), 2 );
  ASSERT_NEAR( xk( 0 ), 0.1, 1.0e-10 );
  ASSERT_NEAR( xk( 1 ), 0.1, 1.0e-10 );
}

TEST( ManifoldMapping, fineModelParabolaEvaluation )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelParabola( 1.0e-7 ) );

  fsi::vector x( 2 ), R( 3 ), output( 3 );
  x << 1, 2;

  fineModel->evaluate( x, output, R );

  ASSERT_NEAR( R( 0 ), 1, 1.0e-10 );
  ASSERT_NEAR( R( 1 ), 1, 1.0e-10 );
  ASSERT_NEAR( R( 2 ), 9, 1.0e-10 );
}

TEST( ManifoldMapping, fineModelSimpleEvaluation )
{
  shared_ptr<SurrogateModel> fineModel( new FineModelSimple( 1.0e-7 ) );

  fsi::vector x( 2 ), R( 3 ), output( 3 );
  x << 1, 2;

  fineModel->evaluate( x, output, R );

  ASSERT_NEAR( R( 0 ), 2, 1.0e-10 );
  ASSERT_NEAR( R( 1 ), 2, 1.0e-10 );
  ASSERT_NEAR( R( 2 ), 4, 1.0e-10 );
}

TEST( ManifoldMapping, coarseModelOptimum )
{
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  fsi::vector x0( 2 ), xk( 2 );
  x0.setZero();

  coarseModel->optimize( x0, xk );

  ASSERT_NEAR( xk( 0 ), 0, 1.0e-10 );
  ASSERT_NEAR( xk( 1 ), 0, 1.0e-10 );
}

TEST( ManifoldMapping, coarseModelOptimum1 )
{
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  fsi::vector y( 3 ), x0( 2 ), xk( 2 );
  x0.setZero();
  y << 0.081, 0.1, 0.121;

  coarseModel->optimize( y, x0, xk );

  ASSERT_NEAR( xk( 0 ), 0.02, 1.0e-2 );
  ASSERT_NEAR( xk( 1 ), 0.101, 1.0e-2 );
}

TEST( ManifoldMapping, coarseModelOptimum2 )
{
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  fsi::vector y( 3 ), x0( 2 ), xk( 2 );
  x0.setZero();
  y << 0.10011, 0.10125, 0.10241;

  coarseModel->optimize( y, x0, xk );

  ASSERT_NEAR( xk( 0 ), 0.001, 1.0e-3 );
  ASSERT_NEAR( xk( 1 ), 0.101, 1.0e-3 );
}

TEST( ManifoldMapping, coarseModelOptimum3 )
{
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  fsi::vector y( 3 ), x0( 2 ), xk( 2 );
  x0.setZero();
  y << 0, -0.4, 0.1;

  coarseModel->optimize( y, x0, xk );

  ASSERT_NEAR( xk( 0 ), 0.05, 1.0e-2 );
  ASSERT_NEAR( xk( 1 ), -0.1, 1.0e-2 );
}

TEST( ManifoldMapping, coarseModelOptimum4 )
{
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  fsi::vector y( 3 ), x0( 2 ), xk( 2 );
  x0.setZero();
  y << 0, -0.35, 0.2;

  coarseModel->optimize( y, x0, xk );

  ASSERT_NEAR( xk( 0 ), 0.1, 1.0e-2 );
  ASSERT_NEAR( xk( 1 ), -0.05, 1.0e-2 );
}

TEST( ManifoldMapping, coarseModelEvaluation )
{
  shared_ptr<SurrogateModel> coarseModel( new CoarseModelParabola() );

  fsi::vector x( 2 ), R( 3 ), output( 3 );
  x << 1, 2;

  coarseModel->evaluate( x, output, R );

  ASSERT_NEAR( R( 0 ), 1, 1.0e-5 );
  ASSERT_NEAR( R( 1 ), 2, 1.0e-5 );
  ASSERT_NEAR( R( 2 ), 3, 1.0e-5 );
}
