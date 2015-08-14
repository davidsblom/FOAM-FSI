#   include "readTimeControls.H"
#   include "readAutoPIMPLEControls.H"
    label oCorr=0;
    lduMatrix::solverPerformance UEqnPerformance;
    lduMatrix::solverPerformance pEqnPerformance;

    bool correctPhi = false;
    if (pimple.found("correctPhi"))
    {
        correctPhi = Switch(pimple.lookup("correctPhi"));
    }

    bool checkMeshCourantNo = false;
    if (pimple.found("checkMeshCourantNo"))
    {
        checkMeshCourantNo = Switch(pimple.lookup("checkMeshCourantNo"));
    }
