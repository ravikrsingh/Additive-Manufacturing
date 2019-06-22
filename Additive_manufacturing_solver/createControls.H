#if defined(NO_CONTROL)
#elif defined(PISO_CONTROL)
    #include "createPisoControl.H"
#elif defined(PIMPLE_CONTROL)
    #include "createPimpleControl.H"
#elif defined(SIMPLE_CONTROL)
    #include "createSimpleControl.H"
#endif


const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");


int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(readScalar(stressControl.lookup("dD")));