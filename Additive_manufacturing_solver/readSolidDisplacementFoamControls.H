nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 10000);
convergenceTolerance = readScalar(stressControl.lookup("D"));