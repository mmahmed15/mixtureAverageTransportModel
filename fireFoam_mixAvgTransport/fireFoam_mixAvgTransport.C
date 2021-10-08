/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    fireFoam

Description
    Transient solver for fires and turbulent diffusion flames with reacting
    particle clouds, surface film and pyrolysis modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiCombustionModel.H"
#include "radiationModel.H"
#include "pimpleControl.H"
#include "SLGThermo.H"
#include "fvOptions.H"

#include "singleStepCombustion.H"
#include "thermoPhysicsTypes.H"
#include "IOmanip.H"

//Mixture Average Transport (MAT)
#include "moleFraction.H"
#include "laminarTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    //#include "infoFieldsOutput.H" // needs revision
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info << "solving two-way coupled soot-gas conservation equations" << endl;
    
        #include "rhoEqn.H"

        //update Mixture Average Transport (MAT)
        transportModel.update();

        // --- PIMPLE loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "YEEqn.H"

            //update mole fractions for MAT
            moleFraction_.update();

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rho = thermo.rho();

        #include "infoOutput.H" // need revision

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"; // kvm

        if (runTime.writeTime()) // kvm
        {
            Info<< " +"; // kvm
        }
        Info<< nl << endl; // kvm

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
