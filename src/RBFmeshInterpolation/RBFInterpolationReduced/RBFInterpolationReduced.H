
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef RBFInterpolationReduced_H
#define RBFInterpolationReduced_H

#include "dictionary.H"
#include "RBFFunctionTG.H"
#include "simpleMatrix.H"
#include "polyMesh.H"
#include "Time.H"

#include "demandDrivenData.H"

#   include "cpuTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

extern "C"
{
// -------------- Factorizing A for solving inv(A) ----------------
// ------- For real symmetric indefinite matrix in packed form ------
void dsptrf_(
    const char * UPLO,
    const int * N,
    const double * A,
    int * ipiv,
    int * info
    );

/*
 * UPLO    (input) CHARACTER*1
 *        = 'U':  Upper triangle of A is stored;
 *        = 'L':  Lower triangle of A is stored.
 *
 * N       (input) INTEGER
 *        The order of the matrix A.  N >= 0.
 *
 * AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
 *        On entry, the upper or lower triangle of the symmetric matrix
 *        A, packed columnwise in a linear array.  The j-th column of A
 *        is stored in the array AP as follows:
 *        if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
 *        if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
 *
 *        On exit, the block diagonal matrix D and the multipliers used
 *        to obtain the factor U or L, stored as a packed triangular
 *        matrix overwriting A (see below for further details).
 *
 * IPIV    (output) INTEGER array, dimension (N)
 *        Details of the interchanges and the block structure of D.
 *        If IPIV(k) > 0, then rows and columns k and IPIV(k) were
 *        interchanged and D(k,k) is a 1-by-1 diagonal block.
 *        If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
 *        columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
 *        is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
 *        IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
 *        interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
 *
 * INFO    (output) INTEGER
 *        = 0: successful exit
 *        < 0: if INFO = -i, the i-th argument had an illegal value
 *        > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
 *             has been completed, but the block diagonal matrix D is
 *             exactly singular, and division by zero will occur if it
 *             is used to solve a system of equations.
 */

// -------------- Solving inv(A) ----------------------------------
// ------- For real symmetric indefinite matrix in packed form ------
void dsptri_(
    const char * UPLO,
    const int * N,
    double * A,
    int * ipiv,
    double * WORK,
    int * info
    );

/*
 * UPLO    (input) CHARACTER*1
 *        Specifies whether the details of the factorization are stored
 *        as an upper or lower triangular matrix.
 *        = 'U':  Upper triangular, form is A = U*D*U**T;
 *        = 'L':  Lower triangular, form is A = L*D*L**T.
 *
 * N       (input) INTEGER
 *        The order of the matrix A.  N >= 0.
 *
 * AP      (input/output) DOUBLE PRECISION array, dimension (N*(N+1)/2)
 *        On entry, the block diagonal matrix D and the multipliers
 *        used to obtain the factor U or L as computed by DSPTRF,
 *        stored as a packed triangular matrix.
 *
 *        On exit, if INFO = 0, the (symmetric) inverse of the original
 *        matrix, stored as a packed triangular matrix. The j-th column
 *        of inv(A) is stored in the array AP as follows:
 *        if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
 *        if UPLO = 'L',
 *           AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
 *
 * IPIV    (input) INTEGER array, dimension (N)
 *        Details of the interchanges and the block structure of D
 *        as determined by DSPTRF.
 *
 * WORK    (workspace) DOUBLE PRECISION array, dimension (N)
 *
 * INFO    (output) INTEGER
 *        = 0: successful exit
 *        < 0: if INFO = -i, the i-th argument had an illegal value
 *        > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
 *             inverse could not be computed.
 */

// -------------- Solving the system A*X = B -> LAPACK------------------------------------
// ----- For symmmetric positive matrices with Cholesky factorization --------------------
void dspsv_(
    const char * UPLO,
    const int * N,
    const int * nrhs,
    double * A,
    int * ipiv,
    double * b,
    const int * ldb,
    int * info
    );

// UPLO : (in) = 'U':  Upper triangle of A is stored;
// = 'L':  Lower triangle of A is stored.
// N    : (in) the number of rows in matrix A (equal to the number of columns)
// nrhs : (in) the number of right hand sides, i.e. the number of columns of b
// A    : (in/out) a (N x (N+1)/2) matrix. On entry, the upper (or lower) triangle of the symmetric
// coefficient matrix A, packed columnwise in a linear array.  The j-th column of A
// is stored in the array AP as follows:
// if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
// if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
// On exit, if INFO = 0, the factor U or L from the Cholesky
// factorization A = U**T*U or A = L*L**T, in the same storage format as A.
// ipiv : (out) vector of lenght N. The pivot indices that define the permutation matrix P;
// row i of the matrix was interchanged with row IPIV(i).
// B    : (in/out) On  entry,  the (N x nrhs) matrix of right hand side matrix B. On exit,
// if INFO = 0, the (N x nrhs) solution matrix X.
// ldb  : (in) The leading  dimension of the array  B, i.e. the number of rows of B.
// info : (out)  = 0:  successful exit
// < 0:  if INFO = -i, the i-th argument had an illegal value
// > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization  has  been completed,
// but the factor U is exactly singular, so the solution could not be  computed.

// ----- computing the matrix-matrix product: C = alpha1*A*B + alpha2*C -> BLAS3 -----------
void dgemm_(
    const char * transa,
    const char * transb,
    const int * M,
    const int * N,
    const int * K,
    const double * alpha1,
    double * A,
    const int * lda,
    double * B,
    const int * ldb,
    const double * alpha2,
    double * C,
    const int * ldc
    );

// transa : (in) = "N", A is used in the computation.
// = "T", AT is used in the computation.
// = "C", AH is used in the computation.
// transb : (in) = "N", B is used in the computation.
// = "T", BT is used in the computation.
// = "C", BH is used in the computation.
// M      : (in) the number of rows in matrix C.
// N      : (in) the number of columns in matrix C.
// K      : (in) If transa="N", the number of columns in matrix A.
// If transa="T" or "C", the number of rows in matrix A, or
// If transa="N", the number of rows in matrix B.
// If transa="T" or "C", the number of columns in matrix B.
// alpha1 : (in) the scaling constant for the matrix-matrix product.
// A      : (in) If transa="N", (M x N) matrix A.
// If transa="T" or "C", the transposed of the (N x M) matrix A is used.
// lda    : (in) the leading dimension of the array specified for A, i.e. the number of rows of A.
// B      : (in) If transa="N", (N x K) matrix B.
// If transa="T" or "C", the transposed of the (K x N) matrix B is used.
// ldb    : (in) the leading dimension of the array specified for B, i.e. the number of rows of B.
// alpha2 : (in) the scaling constant for the matrix.
// C      : (in/out) On entry, the (M x K) matrix C. On exit, the (M x K) matrix C, containing the
// results of the computation.
// ldc    : (in) the leading dimension of the array specified for C, i.e. the number of rows of C.

// ----- computing the matrix-vector product: y = alpha1*A*x + alpha2*y -> BLAS2 -----------
void dgemv_(
    const char * trans,
    const int * M,
    const int * N,
    const double * alpha1,
    double * A,
    const int * lda,
    double * x,
    const int * incx,
    const double * alpha2,
    double * y,
    const int * incy
    );

// trans  : (in) = "N", A is used in the computation.
// = "T", AT is used in the computation.
// = "C", AH is used in the computation.
// M      : (in) the number of rows in matrix C.
// N      : (in) the number of columns in matrix C.
// alpha1 : (in) the scaling constant for the matrix-matrix product.
// A      : (in) If transa="N", (M x N) matrix A.
// If transa="T" or "C", the transposed of the (N x M) matrix A is used.
// lda    : (in) the leading dimension of the array specified for A, i.e. the number of rows of A.
// x      : (in) If transa="N", (N x 1) vector x.
// If transa="T" or "C", the transposed of the (1 x N) vector x is used.
// incx   : (in) the leading dimension of the array specified for x, i.e. the number of rows of x.
// alpha2 : (in) the scaling constant for the matrix.
// y      : (in/out) On entry, the (M x 1) vector y. On exit, the (M x 1) vector y, containing the
// results of the computation.
// incy   : (in) the leading dimension of the array specified for y, i.e. the number of rows of y.

void dgemtx_(
    const int * M,
    const int * N,
    const double * alpha1,
    double * A,
    const int * lda,
    double * x,
    const int * incx,
    double * y,
    const int * incy
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




namespace Foam
{
  // Forward declaration of classes

  /*---------------------------------------------------------------------------*\
  *                       Class RBFInterpolationReduced Declaration
  \*---------------------------------------------------------------------------*/

  class RBFInterpolationReduced
  {
    // Private data

    // - Mesh reference
    const polyMesh & mesh_;

    // - Dictionary
    const dictionary & dict_;

    // - Reference to control points
    const vectorField & controlPoints_;

    // - Rerefence to internal points
    const vectorField & internalPoints_;

    word RBF_;

    // - RBF function
    autoPtr<RBFFunctionTG> RBFfunc_;

    word Dimension_;

    scalar radius_;

    // - Add polynomials to RBF matrix
    bool polyNomials_;

    // - number of reduced columns, needed for size of Hred
    int NRC;

    // - Storage matrix for Hred
    double * Hred;

    // Private Member Functions

    // - Disallow default bitwise copy construct
    RBFInterpolationReduced( const RBFInterpolationReduced & );

    // - Disallow default bitwise assignment
    void operator=( const RBFInterpolationReduced & );

    // - Clear out
    void clearOut();

    void test() const;

public:

    // Constructors

    // - Construct from components
    RBFInterpolationReduced(
      const polyMesh & mesh,
      const dictionary & dict,
      const vectorField & controlPoints,
      const vectorField & internalPoints
      );


    // Destructor

    ~RBFInterpolationReduced();


    // Member Functions

    // - Interpolate for reduced matrix
    template<class Type>
    tmp<Field<Type> > interpolateRed( const Field<Type> & ctrlField ) const;

    // - Move points
    void movePoints();

    // - Create the H matrix such that di=H*d_mov_control
    void createReducedEvaluationMatrix( labelList & movingControlIndex );
  };


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
 #   include "RBFInterpolationReducedTemplates.C"
#endif

#endif

// ************************************************************************* //
