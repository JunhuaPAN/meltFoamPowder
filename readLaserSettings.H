   IOdictionary laserSettings
    (
        IOobject
        (
            "laserSettings",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Reading effective laser absorptivity beta
    dimensionedScalar Efrac(laserSettings.lookup("Efrac"));

    // Reading optical penetration depth
    dimensionedScalar OPD(laserSettings.lookup("OPD"));

    // Reading spot size w
    dimensionedScalar w(laserSettings.lookup("w"));

