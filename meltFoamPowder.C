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
    meltFoamPowder

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

    // Before the time loop begins, read the array of x and y laser coordinates
    // into a variable
    std::ifstream xpoints("xlaser.txt");
    std::ifstream ypoints("ylaser.txt");

    // Store t, x, y in vectors
    std::vector<scalar> tx, ty, x, y;
    scalar tdata, xdata, ydata;
    // Read x laser data
    while(xpoints >> tdata >> xdata)
    {
      tx.push_back(tdata);
      x.push_back(xdata);
    }
    // Read y laser data
    while(ypoints >> tdata >> ydata)
    {
      ty.push_back(tdata);
      y.push_back(ydata);
    }

    // Create scalar for laser intensity
    dimensionedScalar I0 = 2*P/(pi*w*w);

    Info<< "\nStarting time loop\n" << endl;

    // Get maximum z coordinate of mesh
    volScalarField cellz = mesh.C().component(2);

    // Create depth field
    volScalarField depth = max(cellz) - mesh.C().component(2);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Current time value
        scalar tval = runTime.value();

        // Current laser coordinates
        scalar Xlaser = interp(tval,tx,x);
        scalar Ylaser = interp(tval,ty,y);

        // Make distance from laser center a volScalarField
        volScalarField cellx = mesh.C().component(0);
        volScalarField celly = mesh.C().component(1);
        dimensionedScalar X("X",dimensionSet(0,1,0,0,0,0,0),Xlaser);
        dimensionedScalar Y("Y",dimensionSet(0,1,0,0,0,0,0),Ylaser);

        // Square of distance from laser center
        volScalarField R2 = (X-cellx)*(X-cellx) + (Y-celly)*(Y-celly);

        volScalarField laserSource = 1/rho*2.0*P/(pi*w*w)/th*exp(-2.0*R2/(w*w));

        forAll(laserSource,I)
        {
          if(depth[I] > th.value())
            laserSource[I] = 0.0;
        }

        // Create energy depth function
        //volScalarField zfunc = depth/th;
        /*
        forAll(zfunc,I)
        {
          if(depth[I] > th.value())
            zfunc[I] = 0;
        }
        */

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

        // Calculate any extra quantities necessary
        gradT = fvc::grad(T);
        dTdt = fvc::ddt(T);

        // Calculate max/min variable values over simulation duration
        forAll(T, I)
        {
          // Max heating rate
          if(max_dTdt[I] < dTdt[I])
            max_dTdt[I] = dTdt[I];
          // Max cooling rate (min heating rate)
          if(min_dTdt[I] > dTdt[I])
            min_dTdt[I] = dTdt[I];
          // Max temperature gradient
          if(maxMagGradT[I] < mag(gradT[I]))
            maxMagGradT[I] = mag(gradT[I]);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
