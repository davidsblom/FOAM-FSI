#include "fileTranslationalBody.H"
#include "mathematicalConstants.H"
#include "tensor.H"
#include "addToRunTimeSelectionTable.H"
#include "IFstream.H"
#include "IStringStream.H"

using namespace Foam;

defineTypeNameAndDebug(fileTranslationalBody, 0);
addToRunTimeSelectionTable(bodyMotionFunction, fileTranslationalBody, dictionary);


// * * * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //
void fileTranslationalBody::writeState(const scalar& time)
{
	scalar state = getState(time);
	ofBody_ << time << "\t" << state << endl;
}

scalar fileTranslationalBody::getState(const scalar& time)
{
	/*scalar smoothingScalarTrans = 1;
	if(smoothStart_)
	{
		if(time < 1.0/frequency_){
			smoothingScalarTrans = 0.5 - 0.5*Foam::cos(pi*frequency_*time);
			//Info << "smoothingScalarTrans = " << smoothingScalarTrans << ", amplitude = " << smoothingScalarTrans*amplitude_ << endl;
		}
	}

	vector translationVector = smoothingScalarTrans*scaleFactor_*sin(2*pi*frequency_*time)*translationDirection_;*/

	scalar translation = amplitudes_[timeIndex_];

	return translation;
}

tmp<Field<vectorField> > fileTranslationalBody::calculatePosition(const scalar time)
{
	tmp<Field<vectorField> > tPosition(new Field<vectorField>(patchIDs_.size()));
    Field<vectorField>& position = tPosition();

	scalar state = getState(time);
	vector translationVector = state*translationDirection_;

	forAll(patchIDs_,ipatch){
		position[ipatch] = initialPoints_[ipatch];
		position[ipatch] += translationVector;
	}

	return tPosition;
}

void fileTranslationalBody::setPatchIDs(){
	forAll(patchNames_,ipatch){
		label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
		if(patchID < 0){
			FatalErrorIn("void fileTranslationalBody::setPatchIDs()")
		    << "patchName " << patchNames_[ipatch] << " does not exist."
		    << abort(FatalError);
		}else{
			patchIDs_[ipatch] = patchID;
		}
	}
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * //
tmp<Field<vectorField> > fileTranslationalBody::getChildMotion()
{
	//If the localTimeIndex is not equal to mesh time index, update has not been called
	if(mesh_.time().timeIndex() != localTimeIndex())
	{
		WarningIn("fileTranslationalBody::getChildMotion()")
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

void fileTranslationalBody::setInitialPoints(){
	forAll(patchIDs_,ipatch){
		initialPoints_[ipatch] = mesh_.boundaryMesh()[patchIDs_[ipatch]].localPoints();
	}
}

const labelList fileTranslationalBody::getChildMotionPatchIDs() const{
	return patchIDs_;
}

const wordList fileTranslationalBody::getChildMotionPatchNames() const{
	return patchNames_;
}

scalarField fileTranslationalBody::readFile(string fileName) const{
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

	//Check if size is as epected
	scalar Ndata = 1.0/frequency_/deltaT_;
	if(Ndata-translationList.size() > SMALL)
	{
		FatalErrorIn("void fileTranslationalBody::readFile(string fileName)")
		<< "Number of data points read in is not as expected: "
		<< translationList.size() << " != " << Ndata << "Please ensure they are equal!"
		<< abort(FatalError);
	}

	return scalarField(scaleFactor_*translationList);
}

void fileTranslationalBody::checkTimeStepSettings() const
{
	if(mesh_.time().deltaT().value()!=deltaT_)
	{
		FatalErrorIn("void fileTranslationalBody::checkTimeStepSettings()")
		<< "deltaT for calculation (" << mesh_.time().deltaT().value() << ") != deltaT for file ("
		<< deltaT_ << "). Please ensure they are equal!"
		<< abort(FatalError);
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fileTranslationalBody::fileTranslationalBody(const dictionary& dict,const fvMesh& mesh, const word name):
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
deltaT_(readScalar(dict.lookup("deltaT"))),
timeIndex_(0),
smoothStart_(false),
currentPosition_(patchIDs_.size()),
prevPosition_(patchIDs_.size()),
initialPoints_(patchIDs_.size()),
ofBody_("body-"+name+"-state.dat")
{
	checkTimeStepSettings();
	amplitudes_ = readFile(fileName_);

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
		ofBody_ << "time [s]\ttranslation [m]" << endl;
	}
	else
	{
		rm("body-"+name+"-state.dat");
	}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fileTranslationalBody::~fileTranslationalBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Field<vectorField>& fileTranslationalBody::getPosition() const
{
	return currentPosition_;
}

void fileTranslationalBody::update()
{
	if(localTimeIndex() < mesh_.time().timeIndex()){
		forAll(prevPosition_,ipatch){
			prevPosition_[ipatch] = currentPosition_[ipatch];
		}

		timeIndex_++;
		if(timeIndex_>=amplitudes_.size())
		{
			timeIndex_=0;
		}
	}
	bodyMotionFunction::update();//This updates localTimeIndex

	//Since this is a forced motion the position can be calculated here
	currentPosition_ = calculatePosition(mesh_.time().value());
}

void fileTranslationalBody::write()
{
	if(writeToFile())
	{
		writeState(mesh_.time().value());
	}
}
