/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFdD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Reference
    Universität Bayreuth
    Lehrstuhl für Technische Thermodynamik und Trasnportprozesse - LTTT
    Fabian Rösler
    Universitätsstraße 30
    95440 Bayreuth
    Tel.: +49 (921) 55-7163

Application
    tempResidualStressFoam

dDescription
    Solves a convection dominated solid/liquid phase change process.
    Convection is induced by Boussinesq approximation.
    Phase change is described by means of a error-function fit-function.
    The temperature profile is then used to calculate the residual elastic stress

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "pimpleControl.H"

//Added For stress calculation 
#include "Switch.H"
#include "incrementTractionDisplacementFvPatchVectorField.H"
//#include "marangoniFvPatchVectorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //Added For stress calculation 
    #include "postProcess.H"
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    //#include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
   // #include "readTimeControls.H"
    #include "createControls.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setDeltaT.H"
    //#include "marangoniFvPatchVectorField.H"
   // #include "expressionSource.H"
   // pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
   
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"
        
        //--- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
        }
        gradT=fvc::grad(T);
        
        //-- Solve for residual elastic stress
       
            Info<< "\nCalculating displacement field\n" << endl;
            #include "readSolidDisplacementFoamControls.H"
            int iCorr = 0;
            scalar initialResidual = 0;
            #include "solveDisplacementField.H"
            #include "CalculateStress.H"
               runTime.write();
             // Told = T;
        

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
