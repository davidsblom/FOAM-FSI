#include "fourierRotatingBody.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "IStringStream.H"

using namespace Foam;

defineTypeNameAndDebug(fourierRotatingBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, fourierRotatingBody, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void fourierRotatingBody::writeState(const scalar& time)
{
	scalar state = getState(time);
	ofBody_ << time << "\t" << state*180.0/mathematicalConstant::pi << endl;
}

scalar fourierRotatingBody::getState(const scalar& time)
{
	const scalar& pi = mathematicalConstant::pi;

	scalar smoothingScalar = 1;
	if(smoothStart_)
	{
		if(time < 1.0/frequency_){
			smoothingScalar = 0.5 - 0.5*Foam::cos(pi*frequency_*time);
			//Info << "smoothingScalarTrans = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarTrans*amplitude_ << endl;
		}
	}

	scalar state = fourierA0_/2;
	forAll(fourierAn_,iFourier){
		state += 	fourierAn_[iFourier]*cos((iFourier+1)*2.0*pi*frequency_*time) +
						fourierBn_[iFourier]*sin((iFourier+1)*2.0*pi*frequency_*time);
	}
	state *= smoothingScalar;

	return state;
}

tmp<Field<vectorField> > fourierRotatingBody::calculatePosition(const scalar time)
{
	tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

	scalar rotAngle = getState(time);
	tensor RzCur(cos(rotAngle), -sin(rotAngle), 0, sin(rotAngle), cos(rotAngle), 0, 0, 0, 1);
	forAll(patchIDs_,ipatch){
		vectorField rotationField = (RzCur & ( initialPoints_[ipatch] - origin_ )) + origin_;
		position[ipatch] = rotationField;
	}

	return tPosition;
}

void fourierRotatingBody::setPatchIDs(){
	forAll(patchNames_,ipatch){
		label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
		if(patchID < 0){
			FatalErrorIn("void fourierRotatingBody::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > fourierRotatingBody::getChildMotion()
{
	//If the localTimeIndex is not equal to mesh time index, update has not been called
	if(mesh_.time().timeIndex() != localTimeIndex())
	{
		WarningIn("fourierRotatingBody::getChildMotion()")
    		<< "Local time index is not equal to global time index. Update functions needs to be called" << endl;
		update();
	}

	tmp<Field<vectorField> > tMotion(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& motion = tMotion();

    forAll(motion,ipatch){
    	//Set motion
    	motion[ipatch] = currentPosition_[ipatch] - prevPosition_[ipatch];
    	//Update prevPosition
    	prevPosition_[ipatch] = currentPosition_[ipatch];
    }

	return tMotion;
}

void fourierRotatingBody::setInitialPoints(){
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList fourierRotatingBody::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList fourierRotatingBody::getChildMotionPatchNames() const{
	return patchNames_;
}

scalarField fourierRotatingBody::readFile(string fileName) const{
	DynamicList<scalar> amplitudesList(0);

	IFstream fs(fileName);

	Info << "Start finding file " << fs.name() << "....File Found? " << fs.opened() << endl;

	while(true)
	{
		scalar x(0);//scalar for x;
		string lin;
		fs.getLine(lin);//get line from stream;

		if (fs.eof() || lin.length()<3)
		{
			break;
		}

		IStringStream ss(lin);
		ss >> x;//stream coordinates to scalar variables
		amplitudesList.append(x);

	}
	amplitudesList.shrink();

	return scalarField(scaleFactor_*amplitudesList);
}

void fourierRotatingBody::createFourierCoeffs()
{
	const scalar& pi = mathematicalConstant::pi;
	scalar T=1.0/frequency_;
	scalar w=2.0*pi/T;
	scalar dt=T/(amplitudes_.size()-1);
	fourierA0_=2.0/T*(sum(amplitudes_*dt));

	scalarField times(amplitudes_.size(),0);
	forAll(times,itime){
		times[itime] = itime*dt;
	}

	forAll(fourierAn_,iFourier){
		fourierAn_[iFourier]=2.0/T*(sum( (amplitudes_*cos((iFourier+1)*w*times)) * dt));
		fourierBn_[iFourier]=2.0/T*(sum( (amplitudes_*sin((iFourier+1)*w*times)) * dt));
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fourierRotatingBody::fourierRotatingBody(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
fileName_(dict.lookup("filename")),
scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
frequency_(readScalar(dict.lookup("frequency"))),
origin_(dict.lookup("origin")),
inputInDegrees_(readBool(dict.lookup("inputInDegrees"))),
amplitudes_(0),
N_(readLabel(dict.lookup("fourierModes"))),
fourierA0_(0),
fourierAn_(N_,0),
fourierBn_(N_,0),
smoothStart_(false),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
	if(inputInDegrees_)
	{
		scaleFactor_ *= mathematicalConstant::pi/180.0;
	}
	amplitudes_ = readFile(fileName_);
	createFourierCoeffs();

	if(dict.found("smoothStart"))
	{
		smoothStart_ = readBool(dict.lookup("smoothStart"));
	}

	//Set patchIDs_ according to patchNames_
	setPatchIDs();

	//Initialize old and current position
	forAll(patchIDs_,ipatch){
		currentPosition_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
		prevPosition_[ipatch] = currentPosition_[ipatch];
	}

	//Set initial points to mesh points
	setInitialPoints();

	if(writeToFile())
	{
		ofBody_ << "time [s]\trotation [deg]" << endl;
	}
	else
	{
		rm("body-"+name+"-state.dat");
	}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fourierRotatingBody::~fourierRotatingBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& fourierRotatingBody::getPosition() const
{
	return currentPosition_;
}

void fourierRotatingBody::update()
{
	if(localTimeIndex() < mesh_.time().timeIndex()){
		forAll(prevPosition_,ipatch){
			prevPosition_[ipatch] = currentPosition_[ipatch];
		}
	}
	bodyMotionFunction::update();//This updates localTimeIndex

	//Since this is a forced motion the position can be calculated here
	currentPosition_ = calculatePosition(mesh_.time().value());
}

void fourierRotatingBody::write()
{
	if(writeToFile())
	{
		writeState(mesh_.time().value());
	}
}
