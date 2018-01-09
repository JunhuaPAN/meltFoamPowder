/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
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
    meltFoamPowder4a

Description
    Solves a convection dominated solid/liquid phase change process.
    Convection is induced by Boussinesq approximation.
    Phase change is described by means of a error-function fit-function.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mathematicalConstants.H"
#include "pimpleControl.H"

// Use c++ vectors
#include <vector>

// Read file
#include <fstream>

// Custom interpolation function
#include "interp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"  

    // For interpolating a table read OpenFoam-style
    //#include"interpolationTable.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Read the laser path data for solver to process
    #include "readLaserPath.H"

    // Assume initial laser velocity is zero
    scalar Vlaser[2];
    Vlaser[0] = 0.0;
    Vlaser[1] = 0.0;

    // Add code: set tval and deltaT already and calculate xlaserOld from that
    // scalar tval   = runTime.value();
    // scalar deltaT = runTime.deltaTValue();
    // scalar XlaserOld = interp(tval,tx,x);
    // scalar YlaserOld = interp(tval,ty,y);
    // Then, chang the following
    // - Get rid of "scalar" declarations of tval, deltaT, X/YlaserOld

    // Initial laser coordinates
    scalar XlaserOld = x[0];
    scalar YlaserOld = y[0];

    // Declare Xlaser/Ylaser
    scalar Xlaser, Ylaser;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Current time value
        scalar tval = runTime.value();

        // Time step value
        scalar deltaT = runTime.deltaTValue();

        // Current laser coordinates via time-interpolation of input data
        Xlaser = interp(tval,tx,x);
        Ylaser = interp(tval,ty,y);

        // Calculate laser velocity based on old coordinate values
        Vlaser[0] = (Xlaser - XlaserOld)/deltaT;
        Vlaser[1] = (Ylaser - YlaserOld)/deltaT;

        // Check that velocity calculation is accurate
        scalar Vmag = Foam::sqrt(Vlaser[0]*Vlaser[0] + Vlaser[1]*Vlaser[1]);
        Info << "Vx   = " << Vlaser[0] << endl;
        Info << "Vy   = " << Vlaser[1] << endl;
        Info << "Vmag = " << Vmag << nl << endl;

        // Make distance from laser center a volScalarField
        volScalarField cellx = mesh.C().component(0);
        volScalarField celly = mesh.C().component(1);
        dimensionedScalar X("X",dimensionSet(0,1,0,0,0,0,0),Xlaser);
        dimensionedScalar Y("Y",dimensionSet(0,1,0,0,0,0,0),Ylaser);

        // Square of distance from laser center
        volScalarField R2 = (X-cellx)*(X-cellx) + (Y-celly)*(Y-celly);

        // Laser heating source term - turn off after last laser table value
        laserSource = (1.0/rho)*(2.0*P)/(pi*w*w)*edensity*exp(-2.0*R2/(w*w));
        if(tval > tx[tx.size()-1])
          laserSource = 0*laserSource;

        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
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

        // Write max temperature to screen
        Info << "Tmax = " << Foam::max(T) << nl << endl;

        // Calculate any extra quantities necessary
        gradT = fvc::grad(T);
        dTdt = fvc::ddt(T);

        // Generate outputs for solidification quantification
        #include "solidificationCalcs.H"

        // Calculate gradient of alpha and see if it is normal vector
        // to SL interfac
        grad_alpha = fvc::grad(alpha); 

        // Update 'old' laser coordinates for next iteration
        XlaserOld = Xlaser;
        YlaserOld = Ylaser;
	// Deal with things related to alpha and calculating gg etc
        #include "alpha.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
