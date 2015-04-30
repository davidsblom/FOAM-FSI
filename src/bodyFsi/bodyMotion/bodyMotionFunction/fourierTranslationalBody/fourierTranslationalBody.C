#include "fourierTranslationalBody.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "IFstream.H"
#include "IStringStream.H"

using namespace Foam;

defineTypeNameAndDebug(fourierTranslationalBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, fourierTranslationalBody, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void fourierTranslationalBody::writeState(const scalar& time)
{
	vector state = getState(time);
	ofBody_ << time << "\t" << state.z()*180.0/mathematicalConstant::pi << "\t" << state.x() << "\t" << state.y()<< endl;
}

vector fourierTranslationalBody::getState(const scalar& time)
{
	const scalar& pi = mathematicalConstant::pi;

	scalar smoothingScalarTrans = 1;
	scalar smoothingScalarRot = 1;
	if(smoothStart_)
	{
		if(time < 1.0/frequency_){
			smoothingScalarTrans = 0.5 - 0.5*Foam::cos(pi*frequency_*time);
			//Info << "smoothingScalarTrans = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarTrans*amplitude_ << endl;
		}
	}

	vector translationVector =
		smoothingScalarTrans*scaleFactor_*sin(2*pi*frequency_*time)*translationDirection_;

	vector state(translationVector.x(),translationVector.y(),0);

	return state;
}

tmp<Field<vectorField> > fourierTranslationalBody::calculatePosition(const scalar time)
{
	tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

	vector state = getState(time);
	vector translationVector(state.x(),state.y(),0.0);

	forAll(patchIDs_,ipatch){
		position[ipatch] = translationVector;
	}

	return tPosition;
}

void fourierTranslationalBody::setPatchIDs(){
	forAll(patchNames_,ipatch){
		label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
		if(patchID < 0){
			FatalErrorIn("void fourierTranslationalBody::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > fourierTranslationalBody::getChildMotion()
{
	//If the localTimeIndex is not equal to mesh time index, update has not been called
	if(mesh_.time().timeIndex() != localTimeIndex())
	{
		WarningIn("fourierTranslationalBody::getChildMotion()")
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

void fourierTranslationalBody::setInitialPoints(){
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList fourierTranslationalBody::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList fourierTranslationalBody::getChildMotionPatchNames() const{
	return patchNames_;
}

scalarField fourierTranslationalBody::readFile(string fileName) const{
	DynamicList<scalar> translationList(0);

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
		translationList.append(x);

	}
	translationList.shrink();

	return scalarField(scaleFactor_*translationList);
}

void fourierTranslationalBody::createFourierCoeffs()
{
	const scalar& pi = mathematicalConstant::pi;
	scalar T=1.0/frequency_;
	scalar w=2.0*pi/T;
	scalar dt=T/(amplitudes_.size()-1);
	scalar a0=2.0/T*(sum(amplitudes_*dt));

	scalarField times(amplitudes_.size(),0);
	forAll(times,itime){
		times[itime] = itime*dt;
	}

	scalar N=5;
	scalarField an(N,0);
	scalarField bn(N,0);
	forAll(an,i){
	    an[i]=2.0/T*(sum( (amplitudes_*cos((i+1)*w*times)) * dt));
	    bn[i]=2.0/T*(sum( (amplitudes_*sin((i+1)*w*times)) * dt));
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fourierTranslationalBody::fourierTranslationalBody(const dictionary& dict,const fvMesh& mesh, const word name):
bodyMotionFunction(dict,mesh,name),
dict_(dict),
mesh_(mesh),
name_(name),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
fileName_(dict.lookup("filename")),
scaleFactor_(readScalar(dict.lookup("scaleFactor"))),
frequency_(readScalar(dict.lookup("frequency"))),
translationDirection_(dict.lookup("translationDirection")),
amplitudes_(0),
smoothStart_(false),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
	amplitudes_ = readFile(fileName_);
	createFourierCoeffs();

	translationDirection_ /= mag(translationDirection_) + SMALL;
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
		ofBody_ << "time [s]\ttheta [deg]\tx [m]\ty [m]" << endl;
	}
	else
	{
		rm("body-"+name+"-state.dat");
	}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fourierTranslationalBody::~fourierTranslationalBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& fourierTranslationalBody::getPosition() const
{
	return currentPosition_;
}

void fourierTranslationalBody::update()
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

void fourierTranslationalBody::write()
{
	if(writeToFile())
	{
		writeState(mesh_.time().value());
	}
}
