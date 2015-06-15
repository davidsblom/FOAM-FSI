Info<< "Reading solid body properties\n" << endl;

IOdictionary bodyMotionDict
(
    IOobject
    (
        "bodyMotionDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

bodyCollector bodyCol(bodyMotionDict,mesh);

bool forcing = bodyMotionDict.lookupOrDefault("forcing",false);

vector prescribedForces = vector::zero;
scalar prescribedFrequency = 0;
if(forcing){
    Info << "Add forcing to system" << endl;
    prescribedForces = vector(bodyMotionDict.lookup("prescribedForces"));
    prescribedFrequency = readScalar(bodyMotionDict.lookup("prescribedFrequency"));
}
