/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    DPMFoam-extended

Group
    grpLagrangianSolvers

Description
    Transient solver for the coupled transport of a single kinematic particle
    cloud including the effect of the volume fraction of particles on the
    continuous phase.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "DPMIncompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"

#ifdef MPPIC
    #include "basicKinematicCloud.H"
    #define basicKinematicTypeCloud basicKinematicCloud
#else
    #include "basicKinematicCollidingCloud.H"
    #define basicKinematicTypeCloud basicKinematicCollidingCloud
#endif

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for the coupled transport of a"
        " single kinematic particle cloud including the effect"
        " of the volume fraction of particles on the continuous phase."
    );
    argList::addOption
    (
        "cloud",
        "name",
        "specify alternative cloud name. default is 'kinematicCloud'"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        continuousPhaseTransport.correct();
        muc = rhoc*continuousPhaseTransport.nu();

        Info<< "Evolving " << kinematicCloud.name() << endl;
        kinematicCloud.evolve();

        // Update continuous phase volume fraction field
        alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaPhic = alphacf*phic;

        fvVectorMatrix cloudSU(kinematicCloud.SU(Uc));
        volVectorField cloudVolSUSu
        (
            IOobject
            (
                "cloudVolSUSu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(cloudSU.dimensions()/dimVolume, Zero),
            fvPatchFieldBase::zeroGradientType()
        );

        cloudVolSUSu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
        cloudVolSUSu.correctBoundaryConditions();
        cloudSU.source() = Zero;

//         cloudVolSUSu.primitiveFieldRef() =
//             (cloudSU.diag()*Uc() - cloudSU.source())/mesh.V();
//         cloudVolSUSu.correctBoundaryConditions();
//         cloudSU.source() = cloudSU.diag()*Uc();


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UcEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                continuousPhaseTurbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
