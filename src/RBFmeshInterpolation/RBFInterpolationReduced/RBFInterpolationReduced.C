
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "RBFInterpolationReduced.H"
#include "demandDrivenData.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::RBFInterpolationReduced::clearOut()
{}

void Foam::RBFInterpolationReduced::test() const
{
  scalarField weights = RBFfunc_->weights( controlPoints_, controlPoints_[1] );
  scalarField we = RBFfunc_->weights( controlPoints_ );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFInterpolationReduced::RBFInterpolationReduced(
  const polyMesh & mesh,
  const dictionary & dict,
  const vectorField & controlPoints,
  const vectorField & internalPoints
  )
  :
  mesh_( mesh ),
  dict_( dict ),
  controlPoints_( controlPoints ),
  internalPoints_( internalPoints ),
  RBF_( word( dict.lookup( "RBF" ) ) ),
  RBFfunc_( RBFFunctionTG::New( RBF_, dict ) ),
  Dimension_( word( dict.lookup( "Dim" ) ) ),
  radius_( 1 ),
  polyNomials_( readBool( dict.lookup( "polynomials" ) ) ),
  NRC( 0 )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


Foam::RBFInterpolationReduced::~RBFInterpolationReduced()
{
  clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::RBFInterpolationReduced::movePoints()
{
  clearOut();
}

void Foam::RBFInterpolationReduced::createReducedEvaluationMatrix( labelList & movingControlIndex )
{
  int dimtmp = 3;

  if ( Dimension_ == "TwoD" )
  {
    dimtmp = 2;
  }
  else
  if ( Dimension_ == "ThreeD" )
  {
    dimtmp = 3;
  }
  else
  {
    FatalErrorIn( "void Foam::RBFInterpolationReduced::createReducedEvaluationMatrix(labelList& movingControlIndex)" )
    << "Dimension are neither TwoD or ThreeD"
    << abort( FatalError );
  }

  const int dim = dimtmp;       // 3D problem
  int i;
  int j;

  int NC;
  int NColA;

  if ( !polyNomials_ )
  {
    NC = controlPoints_.size();
    NColA = controlPoints_.size();
  }
  else
  if ( polyNomials_ )
  {
    NC = controlPoints_.size() + 1 + dim;
    NColA = controlPoints_.size() + 1 + dim;
  }
  else
  {}

  // ----------- Construct matrix C in packed storage --------------
  // C = | M_11     P_1 |
  // | P_1^T    0   |

  double * C = new double[NC * (NC + 1) / 2];

  for ( i = NC * (NC + 1) / 2; i--; )
  {
    C[i] = 0.0;
  }

  int k = 0;    // integers for row, column and dimension

  // M_11
  // Only upper half
  scalarField weights = RBFfunc_->weights( controlPoints_ );

  for ( i = 0; i < controlPoints_.size(); i++ )
  {
    for ( j = i; j < controlPoints_.size(); j++ )
    {
      C[i + (j + 1) * j / 2] = weights[i + (j + 1) * j / 2];
    }

    // P_1
    if ( polyNomials_ )
    {
      C[i + (controlPoints_.size() + 1) * ( controlPoints_.size() ) / 2] = 1.0;

      for ( k = dim; k--; )
      {
        C[i + (controlPoints_.size() + k + 2) * (controlPoints_.size() + k + 1) / 2] = controlPoints_[i][k];
      }
    }
  }

  // Determine inv(C)
  int INFO = 0;                      // Initialization
  int * IPIV = new int[NC];          // Pivoting
  double * WORK = new double[NC];
  int INFO2 = 0;

  // Facorize the upper matrix and use this in next function
  dsptrf_( "U", &NC, C, IPIV, &INFO );

  if ( INFO != 0 )
  {}

  // Return the inverse of C into C
  dsptri_( "U", &NC, C, IPIV, WORK, &INFO2 );

  if ( INFO2 != 0 )
  {}

  // Construct the reduced inv(C) matrix from previous packed result in FULL form
  NRC = movingControlIndex.size();
  double * CIR = new double[NC * NRC];    // NC rows by NRC columns inverse C matrix

  for ( j = 0; j < NRC; j++ ) // loop over all reduced columns
  {
    for ( i = 0; i < NC; i++ )      // loop over all rows
    {
      int q = movingControlIndex[j];

      // Use symmetric definition c_ij=c_ji
      int ii = min( i, q );
      int jj = max( i, q );

      CIR[i + j * NC] = C[ii + (jj + 1) * jj / 2];       // CHECK - SIGN!!!!!
    }
  }

  delete[] C;
  delete[] IPIV;
  delete[] WORK;

  // Build up complete Hred matrix row by row
  Hred = new double[internalPoints_.size() * NRC];

  for ( i = 0; i < internalPoints_.size(); i++ )
  {
    double * A = new double[NColA];

    for ( j = NColA; j--; )
    {
      A[j] = 0.0;
    }

    weights = RBFfunc_->weights( controlPoints_, internalPoints_[i] );
    j = 0;
    k = 0;

    // M_21
    for ( j = 0; j < controlPoints_.size(); j++ )
    {
      A[j] = weights[j];
    }

    // P_2
    if ( polyNomials_ )
    {
      A[controlPoints_.size()] = 1;

      for ( k = 0; k < dim; k++ )
      {
        A[controlPoints_.size() + k + 1] = internalPoints_[i][k];
      }
    }

    // coefficients for solving matrix-vector product with BLAS2
    const double alpha1 = 1;
    const double alpha2 = 0;
    int incx = 1;
    int incy = 1;

    // Constructing row of H matrix
    double * Ht = new double[NRC];

    for ( j = 0; j < NRC; j++ )
    {
      Ht[j] = 0.0;
    }

    dgemv_( "T", &NC, &NRC, &alpha1, CIR, &NC, A, &incx, &alpha2, Ht, &incy );

    for ( j = 0; j < NRC; j++ )
    {
      Hred[i + j * internalPoints_.size()] = Ht[j];
    }

    // End constructing row of H matrix and saved into storage

    delete[] A;
    delete[] Ht;
  }

  delete[] CIR;
}

// ************************************************************************* //
