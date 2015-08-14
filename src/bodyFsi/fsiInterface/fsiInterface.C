/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

edit by HB 20090925: adjust for different boundary set.

\*---------------------------------------------------------------------------*/

#include "fsiInterface.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

void fsiInterface::setRelexationParamter()
{

	if(couplingScheme_=="Aitken")
	{
		if(fsiIter_>1)
		{
			relaxFsi_ = -relaxFsiPrev_
				* (gSum(residualPrev_ & (residual_ - residualPrev_)))
				/ (gSum((residual_ - residualPrev_) & (residual_ - residualPrev_)));

			relaxFsi_ = min(relaxFsi_,relaxFsiUpper_);
			relaxFsi_ = max(relaxFsi_,relaxFsiLower_);
		}
		Info << "Fsi Relaxation set to " << relaxFsi_ << " with Aitken" << endl;
	}
	else
	{
		relaxFsi_ = relaxFsi0_;
		Info << "Fsi Relaxation set to " << relaxFsi_ << endl;
	}
}

void fsiInterface::setInterfaceForce(const dimensionedScalar& rhoFluid, const dimensionedScalar& nu,const volVectorField& U,const volScalarField& p)
{
	// Set fluidForces
	fluidForcesPrev_ = fluidForces_;

	fluidPatchForces_ = forcesObj_.calcForceFields();
	label counter = 0;
	forAll(fluidPatchForces_,ipatch){
		forAll(fluidPatchForces_[ipatch],iface){
			fluidForces_[counter] = fluidPatchForces_[ipatch][iface];
			counter++;
		}
	}
}

void fsiInterface::syncPatchForcesAfterRelax()
{
	label counter=0;
	forAll(fluidPatchForces_,ipatch){
		forAll(fluidPatchForces_[ipatch],ivalue){
			fluidPatchForces_[ipatch][ivalue] = fluidForces_[counter];
			counter++;
		}
	}
}

void fsiInterface::relax()
{
	// Set fluidForces with under relaxation
	setRelexationParamter();
	fluidForces_ = fluidForcesPrev_ + relaxFsi_*residual_;// F(k+1) = w*F~(k+1) + (1-w)*F(k) = F(k) + w*(F~(k+1)-F(k)) = F(k) + w*res
	syncPatchForcesAfterRelax();
}

scalar fsiInterface::normResidual()
{
	//relative residual: ||f(k) - f(k-1)||/||f(k)||
	scalar normRes = (::sqrt(gSum(magSqr(residual_))))/(::sqrt(gSum(magSqr(fluidForces_))));

	/*scalar normRes=0;
	label totalPatchSize = residual_.size();
	reduce(totalPatchSize,sumOp<label>());
	normRes = (::sqrt(gSum(magSqr(residual_))))/normValue_/totalPatchSize;*/

	return normRes;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

fsiInterface::fsiInterface
(
	const fvMesh& mesh,
	const IOdictionary& dict,
	const labelList& movingPatchIDs
):
dict_(dict),
mesh_(mesh),
fsiIter_(0),
minFsiIter_(readLabel(dict_.lookup("minIter"))),
maxFsiIter_(readLabel(dict_.lookup("maxIter"))),
minResidual_(readScalar(dict_.lookup("minResidual"))),
relaxFsiUpper_(1.0),
relaxFsiLower_(0.0),
relaxFsi0_(readScalar(dict_.lookup("fsiRelaxation"))),
relaxFsi_(relaxFsi0_),
relaxFsiPrev_(relaxFsi0_),
couplingScheme_(dict_.lookup("couplingScheme")),
writeToFile_(dict_.lookupOrDefault("writeToFile",false)),
reuseRelaxFsi_(dict_.lookupOrDefault("reuseRelaxFsi",false)),
fluidPatchForces_(mesh.boundaryMesh().size(),vectorField(0)),
fluidForces_(0,vector::zero),
fluidForcesPrev_(0,vector::zero),
residual_(0,vector::zero),
residualPrev_(0,vector::zero),
normRes_(0),
ofFsiProp_("fsiProperties.dat"),
movingPatchIDs_(movingPatchIDs),
forcesObj_(word("forces"),mesh,dict.subDict("forces"),movingPatchIDs_)
{
	Info << "FsiInterface uses following patches as interface patches: " << movingPatchIDs_ << endl;

    //- Create initial fluidForces
    label fluidForceSize=0;
    forAll(movingPatchIDs_,ipatch){
    	fluidForceSize += mesh_.boundaryMesh()[movingPatchIDs_[ipatch]].size();
    }
    fluidForces_ = vectorField(fluidForceSize,vector::zero);

    //- Check if upperbound is specified. Otherwise = 1.0
    if(dict_.found("fsiRelaxationUpper")){
    	relaxFsiUpper_ = readScalar(dict_.lookup("fsiRelaxationUpper"));
    	Info << "Upper bound for Aitkens relaxation paramter set to " << relaxFsiUpper_ << endl;
    }
	//- Check if lowerbound is specified. Otherwise = 0.0
    if(dict_.found("fsiRelaxationLower")){
    	relaxFsiLower_ = readScalar(dict_.lookup("fsiRelaxationLower"));
    	Info << "Lower bound for Aitkens relaxation paramter set to " << relaxFsiLower_ << endl;
    }

	//Create header in output file or delete the file
	if(writeToFile_)
	{
		ofFsiProp_ << "time\tnIters\tresidual" << endl;
	}
	else
	{
		rm("fsiProperties.dat");
	}
}

fsiInterface::~fsiInterface()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void fsiInterface::initializeForce(const dimensionedScalar& rhoFluid, const dimensionedScalar& nu,const volVectorField& U,const volScalarField& p)
{
	//From constructor this should be zero size, otherwise is set via setInterfaceForce
	if(fluidForcesPrev_.size()>0)
	{
		WarningIn("fsiInterface::initializeForce()")
			<< "Force already set!";
	}
	else
	{
		//- Set forces
		setInterfaceForce(rhoFluid,nu,U,p);
	}
}

const Field<vectorField>& fsiInterface::getPatchForces() const
{
	//Warning in parallel this only for current processor
	return fluidPatchForces_;
}

const vectorField& fsiInterface::getForces() const
{
	//Warning in parallel this only for current processor
	return fluidForces_;
}

vector fsiInterface::getTotalForce() const
{
	return gSum(fluidForces_);
}

//- Reset variables when new time is entered
void fsiInterface::reset()
{
	fsiIterTotal_ += fsiIter_;
	fsiIter_ = 0;
	if(!reuseRelaxFsi_){
		relaxFsi_ = relaxFsi0_;
		relaxFsiPrev_ = relaxFsi0_;
	}

	residual_ = vectorField(0,vector::zero);
	residualPrev_ = vectorField(0,vector::zero);
}

void fsiInterface::update(const dimensionedScalar& rhoFluid, const dimensionedScalar& nu,const volVectorField& U,const volScalarField& p)
{
	//- Set forces
	setInterfaceForce(rhoFluid,nu,U,p);

	// Set residual according to new (non relaxated) forces
	residualPrev_ = residual_;
	residual_ = fluidForces_ - fluidForcesPrev_;

	// Set prev relaxationParameter to the current one
	relaxFsiPrev_ = relaxFsi_;

	// Increase iter
	fsiIter_++;
}

bool fsiInterface::converged()
{
	normRes_ = normResidual();
	Info << "fsi SubIteration " << fsiIter_ << " ("<<minFsiIter_<<"<->"<<maxFsiIter_ << "): " << normRes_ << "(" << minResidual_ << ")" << nl << endl;
	if((fsiIter_>=maxFsiIter_ || normRes_ < minResidual_) && fsiIter_ >= minFsiIter_)
	{
		return true;
	}
	else
	{
		relax();//relaxes the force for next subiteration
		return false;
	}
}

label fsiInterface::iter() const
{
	return fsiIter_;
}

void fsiInterface::write()
{
	if(writeToFile_)
	{
		ofFsiProp_ << mesh_.time().value() << "\t" << fsiIter_ << "\t" << fsiIterTotal_ << "\t" << normRes_ << endl;
	}
}

} // End namespace Foam

// ************************************************************************* //
