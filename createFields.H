    // Reading fields
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create depth field as (max(z) - z) on the mesh (CGS units)
    volScalarField cellz(mesh.C().component(2));
    volScalarField depth(max(cellz) - cellz);

    // Create the field gradT
    volVectorField gradT
    (
        IOobject
        (
            "gradT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(T) 
    );

    // Create the field partial(T)/partial(t)
    volScalarField dTdt
    (
        IOobject
        (
            "dTdt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(T) 
    );

    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha_vap\n" << endl;
    volScalarField alpha_vap
    (
        IOobject
        (
            "alpha_vap",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Gradient of alpha - for finding normal vector
    volVectorField grad_alpha
    (
        IOobject
        (
            "grad_alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(alpha) 
    );


    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    // Reading transport properties
    Info<< "Reading thermophysical properties\n" << endl;
    #include "readTransportProperties.H"

    // Read Laser settings for this run
    Info<< "Reading laser settings\n" << endl;
    #include "readLaserSettings.H"

    // Calculating fit-parameters for phase change function
    Info<< "Calculating phase change properties\n" << endl;
    
    dimensionedScalar Tmelt
    (
        "Tmelt", 
        (Tl+Ts)/2.0
    );
    
    dimensionedScalar Tvap
    (
        "Tvap", 
        (Tv2+Tv1)/2.0
    );

    // For viscosity calcs, T should be restricted to T > Tl
    volScalarField Tvisc
    (
        IOobject
        (
            "Tvisc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Foam::max(Tl,T)
    );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - Foam::max(beta*(T - Tl), 0.0)
    );

    // D'arcy-type source term field
    volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh
        ),
        DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs)
    );

    // Heat capacity field
    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh
        ),
        b0 + b1*T
    );

    // Kinematic viscosity field
    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh
        ),
        mu0*exp(A/Tvisc)/rho
    );

    // Add dynamic field so it can be accessed by Marangoni BC
    volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh
        ),
        nu*rho
    );

    // Add powder flag variable
    // 1 = powder, 0 = solid
    volScalarField pflag
    (
        IOobject
        (
            "pflag",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Laser heat source field
    volScalarField edensity
    (
        IOobject
        (
            "edensity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Efrac/OPD*exp(-depth/OPD)
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

// Read grain fraction variable gg
Info<< "Reading field gg\n" << endl;
volScalarField gg
(
    IOobject
    (
        "gg",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read radius of the grains
Info<< "Reading field rad\n" << endl;
volScalarField rad
(
    IOobject
    (
        "rad",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read texture of the grains
Info<< "Reading field theta\n" << endl;
volScalarField theta
(
    IOobject
    (
        "theta",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Temperature gradient during solidification
// Units: [T]/[L]
Info<< "Reading field Gn" << endl;
volScalarField Gn
(
    IOobject
    (
        "Gn",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Cooling rate during solidification
// Units: [T]/[time]
Info<< "Reading field CR" << endl;
volScalarField CR
(
    IOobject
    (
        "CR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Solidification rate
// Units: [L]/[time}
Info<< "Reading field Rdot" << endl;
volScalarField Rdot
(
    IOobject
    (
        "Rdot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Gn/Rdot
// Units: ([T]/[L])/([L]/[time]) = [T][time]/[L^2]
Info<< "Reading field Sratio" << endl;
volScalarField Sratio 
(
    IOobject
    (
        "Sratio",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
