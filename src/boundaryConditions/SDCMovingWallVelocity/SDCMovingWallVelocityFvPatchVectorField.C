
#include "SDCMovingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "ClenshawCurtis.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    SDCMovingWallVelocityFvPatchVectorField::SDCMovingWallVelocityFvPatchVectorField(
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF
        )
        :
        fixedValueFvPatchVectorField( p, iF ),
        k( 0 ),
        corrector( false )
    {
        init();
    }

    SDCMovingWallVelocityFvPatchVectorField::SDCMovingWallVelocityFvPatchVectorField(
        const SDCMovingWallVelocityFvPatchVectorField & ptf,
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF,
        const fvPatchFieldMapper & mapper
        )
        :
        fixedValueFvPatchVectorField( ptf, p, iF, mapper ),
        k( 0 ),
        corrector( false )
    {
        init();
    }

    SDCMovingWallVelocityFvPatchVectorField::SDCMovingWallVelocityFvPatchVectorField(
        const fvPatch & p,
        const DimensionedField<vector, volMesh> & iF,
        const dictionary & dict
        )
        :
        fixedValueFvPatchVectorField( p, iF ),
        k( 0 ),
        corrector( false )
    {
        fvPatchVectorField::operator=( vectorField( "value", dict, p.size() ) );

        init();
    }

    SDCMovingWallVelocityFvPatchVectorField::SDCMovingWallVelocityFvPatchVectorField( const SDCMovingWallVelocityFvPatchVectorField & pivpvf )
        :
        fixedValueFvPatchVectorField( pivpvf ),
        k( 0 ),
        corrector( false )
    {
        init();
    }

    SDCMovingWallVelocityFvPatchVectorField::SDCMovingWallVelocityFvPatchVectorField(
        const SDCMovingWallVelocityFvPatchVectorField & pivpvf,
        const DimensionedField<vector, volMesh> & iF
        )
        :
        fixedValueFvPatchVectorField( pivpvf, iF ),
        k( 0 ),
        corrector( false )
    {
        init();
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void SDCMovingWallVelocityFvPatchVectorField::updateCoeffs()
    {
        if ( updated() )
        {
            return;
        }

        assert( timeIntegrationScheme );

        const fvMesh & mesh = dimensionedInternalField().mesh();

        const fvPatch & p = patch();
        const polyPatch & pp = p.patch();

        // Get wall-parallel mesh motion velocity from geometry
        tmp<vectorField> Up =
            pp.faceCentres() / mesh.time().deltaT().value();

        int dim = mesh.nGeometricD();

        fsi::vector rhs( Up->size() * dim ), f( Up->size() * dim ), result( Up->size() * dim ), qold;

        if ( not corrector && k == 0 )
        {
            const pointField & oldPoints = mesh.oldPoints();

            vectorField oldFc( pp.size() );

            forAll( oldFc, i )
            {
                oldFc[i] = pp[i].centre( oldPoints );

                for ( int j = 0; j < dim; j++ )
                    result( i * dim + j ) = oldFc[i][j];
            }
            timeIntegrationScheme->setOldSolution( mesh.time().timeIndex(), result );
        }

        timeIntegrationScheme->getSourceTerm( corrector, k, sweep, mesh.time().deltaT().value(), rhs, qold );
        forAll( Up(), i )
        {
            for ( int j = 0; j < dim; j++ )
            {
                Up()[i][j] -= ( qold( i * dim + j ) + rhs( i * dim + j ) ) / mesh.time().deltaT().value();
                f( i * dim + j ) = Up()[i][j];
                result( i * dim + j ) = pp.faceCentres()[i][j];
            }

            if ( dim == 2 )
                Up()[i][2] = 0;
        }
        timeIntegrationScheme->setFunction( k, f, result );

        const volVectorField & U =
            mesh.lookupObject<volVectorField>
            (
            dimensionedInternalField().name()
            );

        scalarField phip =
            p.patchField<surfaceScalarField, scalar>( fvc::meshPhi( U ) );

        tmp<vectorField> n = p.nf();
        const scalarField & magSf = p.magSf();
        tmp<scalarField> Un = phip / (magSf + VSMALL);

        // Adjust for surface-normal mesh motion flux
        vectorField::operator=( Up() + n() *( Un() - ( n() & Up() ) ) );

        fixedValueFvPatchVectorField::updateCoeffs();
    }

    void SDCMovingWallVelocityFvPatchVectorField::write( Ostream & os ) const
    {
        fvPatchVectorField::write( os );
        writeEntry( "value", os );
    }

    void SDCMovingWallVelocityFvPatchVectorField::init()
    {
        string filename = static_cast<std::string>( this->db().time().rootPath() ) + "/" + static_cast<std::string>( this->db().time().caseName() ) + "/constant/fsi.yaml";

        if ( Pstream::nProcs() > 1 )
            filename = static_cast<std::string>( this->db().time().rootPath() ) + "/" + static_cast<std::string>( this->db().time().caseName() ) + "/../constant/fsi.yaml";

        YAML::Node config = YAML::LoadFile( filename );

        assert( config["sdc"] || config["esdirk"] );

        if ( config["sdc"] )
        {
            YAML::Node sdcConfig( config["sdc"] );
            assert( sdcConfig["convergence-tolerance"] );
            assert( sdcConfig["number-of-points"] );
            assert( sdcConfig["quadrature-rule"] );

            int n = sdcConfig["number-of-points"].as<int>();
            double tol = sdcConfig["convergence-tolerance"].as<double>();
            std::string quadratureRule = sdcConfig["quadrature-rule"].as<std::string>();

            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;

            if ( quadratureRule == "gauss-radau" )
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussRadau<scalar>( n ) );

            if ( quadratureRule == "gauss-lobatto" )
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::GaussLobatto<scalar>( n ) );

            if ( quadratureRule == "clenshaw-curtis" )
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::ClenshawCurtis<scalar>( n ) );

            if ( quadratureRule == "uniform" )
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >( new fsi::quadrature::Uniform<scalar>( n ) );

            timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::SDC( quadrature, tol ) );
        }

        if ( config["esdirk"] )
        {
            YAML::Node esdirkConfig( config["esdirk"] );

            assert( esdirkConfig["method"] );

            std::string method = esdirkConfig["method"].as<std::string>();

            timeIntegrationScheme = std::shared_ptr<sdc::TimeIntegrationScheme> ( new sdc::ESDIRK( method ) );
        }
    }

    void SDCMovingWallVelocityFvPatchVectorField::setSDCInfo(
        bool corrector,
        int k,
        int sweep
        )
    {
        this->k = k;
        this->sweep = sweep;
        this->corrector = corrector;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField
    (
        fvPatchVectorField,
        SDCMovingWallVelocityFvPatchVectorField
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// ************************************************************************* //
